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
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt
size=${datadir}/genome/sizes.genome #size of chrm
##
cd $dir

#source miniconda3/bin/activate
conda activate basic_tools

for SAMPLE_ID in `cat sample.txt`; do
##sort
#sorting is really important and you will do it many times in the same pipeline usually. A lot of the programs
    #we use expect files to be sorted and break when they are not sorted. Thefore, depending on the file type,
    #bedtools sort or samtools sort for bed files and sam/bam files, respectively, will be very handy.
  samtools sort ${SAMPLE_ID}_chromap.bam -@ ${NSLOTS}  -o ${SAMPLE_ID}_chromap_sorted.bam
#the -@ flag tells samtools how many threads are available for parallelization.
    #many but not all programs will be able to do this (famously, bedtools does not)
    #the only way to know, and to know which flag to use to provide them with the number of cores,
    #is to google it.
    #the ${NSLOTS} variable is AGE specific, and it is a variable that is automatically created when you submit a job
    #and holds the number of processors available to you based on what you requested.

#filtering out mitochondrial reads
## Hi don't run me.

#samtools index ${sorted.bam.file} 

#the above command makes an index file for the bam file, it will be called the same name but .bai will be appended to it.

#samtools idxstats ${sorted.bam.file} | cut -f1 | grep -v Mt | xargs samtools view -b ${sorted.bam.file}  > ${sorted-noMT.bam.file}

#the above command does a few things piped along so we'll split it by pipe.
#first command prints out statistics for the bam file based on the file and the inex file created previously.
#the second command takes that output and cuts just the firsty entry per line, which will be the chromosome name.
#the third command is a grep call, which is a general "search for this text and print the line" unix utility.
    #it is searching for "Mt" which is the name of the mitochondrial chromosome, and the -v flag tells it to print
    #everything that does NOT match that text. If you want to filter more than one chromosome just stack the pipes, like
    # grep -v Mt | grep -v chrX | grep -v chrY etc etc
#xargs is another unix utility that essentially submits a set of commands with one argument at a time and handles the results on the back end.
    #basically what it is doing is taking the piped results, which will be a list of the headings you did not remove, in this case everything but
    #the mitochondria, and submits them to samtools view, then prints them to a new bam file one heading at a time.


##sort
  samtools sort ${SAMPLE_ID}_chromap_sorted.bam -@ ${NSLOTS}  -o ${SAMPLE_ID}_treat.bam
##convert to bw
  samtools index ${SAMPLE_ID}_treat.bam
  conda activate deeptools
  bamCoverage -p max -b ${SAMPLE_ID}_treat.bam  --normalizeUsing RPKM  -v  -o ${SAMPLE_ID}_norm.bw ##you can use this on the genome browser
  #normalizing is usually a good idea so you can compare the bigwigs to each other, RPKM is the standard method people use but there are others
    #with strengths and weaknesses, google will help you here.
  conda activate basic_tools
done