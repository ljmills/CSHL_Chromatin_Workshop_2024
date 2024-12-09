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

#to convert to bam you will need to know the sizes of each chromosome.

datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
cd ${datadir}/genome
conda activate basic_tools

#samtools faidx will make an index for the genome's fasta file.
samtools faidx hg38_chr22.fa
#cut is a unix utility that takes just the first and second fields in each line (in this case)
    #and outputs them to the sizes.genome file, which will be formatted:
        #chr1   12345
        #chr2   12345
        #chr3   12345
cut -f1,2 hg38_chr22.fa.fai > sizes.genome
#this just prints it to make sure it worked well.
cat sizes.genome