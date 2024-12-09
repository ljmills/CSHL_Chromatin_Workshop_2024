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

#source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate macs2

input_combinations_broad=(
"SRR5063143_naive_H3K27ac_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063144_naive_H3K27ac_treat.bam,SRR5063154_naive_input_treat.bam"
)

input_combinations_narrow=(
"SRR5063149_naive_H3K4me3_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063150_naive_H3K4me3_treat.bam,SRR5063154_naive_input_treat.bam"
)


for files in ${input_combinations_broad[@]}; do

#the below line is handy for reading in tables, for example. In this case the "table" is generated
#in the previous lines, but you could imagine them written as a comma-separated file like a .csv
    #the <<< is basically telling unix to run whatever is to the left of it on whatever is to the right of it.
    #IFS stands for Internal Field Separator, and is just a way to tell unix the delimiter it should use to split
    #the input. read -r is used to take input from the user and assign it to a variable (or more).
    #because you are running it before the <<< you are telling read to treat $files as the user input and assign it to
    #variables $file1 and $file2, and to split it based on the IFS of ','

IFS=',' read -r file1 file2 <<< $files
##macs2-broad
macs2 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g hs --nomodel --shift -100 --extsize 200 \
        -n ${file1%_treat.bam} --broad \
        --outdir . 2> ${file1%_treat.bam}_broad_macs2.log
done

#macs2 will call peaks based on your bam(s). for more detail you can look at the macs manuals.
    #in this specific case we are telling it to use the human genome size with -g hs, we are telling it
    #not to make a model of the binding, because there are some assays that result in reads offset from
    #the actual site of binding for whatever you are testing for, but chip seq on histone modifications
    #is not one of them. (chip seq on transcription factors is, though)
    #instead we are going to manually shift and extend the read to make up for the fact that we are using single end
    #data instead of paired end data, so that we can have a longer insert sizes around the place we are sequencing.
    #--outdir . is just a way unix uses to say "the current directory" or in the case of moving files "the current file name."
    #the 2> *.log is telling the command to redirect the error logs to that file.

    #In this specific case, we need to provide the -c control bam, but not all assays require this. All chip seq does but
    #for example ATACseq does not. If you simply do no include a -c flag and file it will not use one.

for files in ${input_combinations_narrow[@]}; do

IFS=',' read -r file1 file2 <<< $files
##macs2-narrow
macs2 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g hs  --shift -100 --extsize 200 \
        -n ${file1%_treat.bam}  \
        --outdir . 2> ${file1%_treat.bam}_narrow_macs2.log

done

#MACS outputs the peak files in .narrowPeak or .broadPeak format.
#I would usually recommend clearly keeping the labels on these so you don't lose track,
#but for downstream convenience we will rename them today.
for peakfile in `ls *.broadPeak`
do
mv ${peakfile} ${peakfile//\.broadPeak/\.Peak}
done

for peakfile in `ls *.narrowPeak`
do
mv ${peakfile} ${peakfile//\.narrowPeak/\.Peak}
done

