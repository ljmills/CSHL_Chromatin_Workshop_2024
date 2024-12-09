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

#THESE ARE INCOMPLETE PATHS, YOU NEED TO ADD TO THEIR BEGINNING
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir=${datadir}/data/subset
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt
size=${datadir}/genome/sizes.genome #size of chrm
##

cd $dir

# source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate basic_tools

#list was previously set in the variables at the beginning as "sample.txt" and is a file that contains library names one per line.
# the backticks tells unix to run whatever is inside them and replace the spot in the line with the output.
# it then assings those to SAMPLE_ID, once for each time in the loop.
#for each of those it uses bedtools to convert them to bam.
for SAMPLE_ID in `cat $list`; do
#convert alignments to BAM
        bedtools bedtobam  -i ${SAMPLE_ID}_chromap.bed -g $size > ${SAMPLE_ID}_chromap.bam #input files are the same from chromap
done