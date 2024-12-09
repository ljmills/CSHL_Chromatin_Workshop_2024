#!/bin/bash
#$ -S /bin/bash
#$ -pe threads 1
#$ -l mfree=120G -l h_rt=100:00:00
################################################################
##select random fastq - subsample
##############################################################

##
out="subset"
##
mkdir -p $out

# source ~/miniconda3/bin/activate
conda activate sra_tools

for i in *.fastq; do
echo $i
seqtk sample -s100 $i 3000000 > $out/${i}.fastq
done

#this will "randomly" take 3 million reads from each file that matches *.fastq and deposit it
#in the $out directory with the same name.

#a note: the custom names that are more informative present in the later scripts were added by hand using the mv command.