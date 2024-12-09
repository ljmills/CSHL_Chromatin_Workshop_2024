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

#The things above this line are an AGE header for submitting the job with qsub.

#THESE ARE NOT ABSOLUTE PATHS YOU NEED TO ADD TO THEM
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir="${datadir}/data/subset"
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt

#you will need to uncomment and run this line to launch conda, so that you can load the conda environments
#source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate chromap

cd $dir
for SAMPLE_ID in `cat $list`; do
# map using chromap, output is bed file
chromap --preset chip -x $index -r $genome -q 20 --min-read-length 10   -1  ${SAMPLE_ID}.fastq  -o  ${SAMPLE_ID}_chromap.bed
done
conda activate basic_tools
#dos2unix basically just converts all the End Of Line (EOL) marker from DOS/Windows format to unix, so from \r\n to \n
#the reason this is here is just because when I tried to run "convert to bam" later, there was an error that was fixed by this.
dos2unix *.bed