#!/bin/bash
#$ -l mfree=5G
#$ -pe threads 8
#$ -l h_rt=100:00:00
#$ -cwd

#what the above commands are doing:
#first line is shebang, it tells the environment to read this as a shell script
#second line requests resources (-l) of 5G PER CORE
#third line requests the processing environment (-pe) to be 8 threads.
    #this will vary between institution. At my home institution, for example, the syntax is '-pe serial 8'
    #memory is counted per core, so in this case the job is requesting 40G total for use by 8 cores.
#the fourth line is asking to be allowed to use 100 hours for this job. Depending on your institution you will have
    #different requirements for declaring the time it will run. If the job finshes early it will terminate, so this
    #actually just means a time at which the job will be unceremoniously killed if it is still running.
#the fifth line is saying to run in the "current working environment" which is handy to include because it means the
#error logs will automatically be output into the directory you are in, and it will also load whatever gets loaded on
#your directory when you launch it (like if you have conda launch on startup for example).

#####################################################################################################################
#The things above this line are an AGE header for submitting the job with qsub.

datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir=${datadir}/data/subset

cd ${dir}

mkdir -p consensus_matrixes/
mkdir -p deeptools_graphs/

#source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate deeptools

input_combinations=(
"SRR5063143_naive_H3K27ac_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063144_naive_H3K27ac_treat.bam,SRR5063154_naive_input_treat.bam"
"SRR5063149_naive_H3K4me3_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063150_naive_H3K4me3_treat.bam,SRR5063154_naive_input_treat.bam"
)

for files in ${input_combinations[@]}; do

IFS=',' read -r file1 file2 <<< $files
#for more information on the above line check the call_peaks.sh script.
bamCompare -b1 ${file1} -b2 ${file2} -o ${file1//_treat\.bam/_differential}.bw
#bamCompare is a deeptools command that makes a bigwig by essentially subtracting one bam
    #from the other.
peakfile=${file1//_treat.bam/_peaks\.Peak}
bw_file=$(echo ${file1//_treat\.bam/_differential}.bw)
#these just set vairables to point to relevant files for ease and readability

computeMatrix scale-regions -p ${NSLOTS} -S ${bw_file} -R ${peakfile} -b 3000 -a 3000 \
        -o consensus_matrixes/${peakfile//_peaks\.Peak/\.matrix}
#this takes one or more bigwig files and one or more bed files, and essentially stakcs and scales the bed file regions horizontally
    #so they all appear the same size
    #so they are all the same width, and calculates the depth at a bin sliding across the bed region
    #the -b 3000 command tells computeMatrix to expand the bed regions upstream by 3000 bp,
    #the -a does the same thing but downstream of the bed region.
    #-o determines the output as usual.
    #one note: this script is VERY slow if not paralellized, so make sure to include a -p ${NSLOTS} argument.

plotHeatmap -m consensus_matrixes/${peakfile//_peaks\.Peak/\.matrix} -o deeptools_graphs/${peakfile//_peaks\.Peak/\.pdf} \
        --dpi 300 --startLabel "Peak Start" --endLabel "Peak End" -x "Distance" --heatmapWidth 12 --regionsLabel "Peaks"

#this guy makes the heatmap from the consensus matrix that was just made in the previous line. --dpi is dots per inch, it
    #essentially controls the quality. the others are all labels and such for customizing the figure.
    #the suffix on the -o will determine the format for output.

done
