#!/bin/bash
#$ -S /bin/bash
#$ -pe threads 1
#$ -l m_mem_free=80G -l h_rt=170:00:00
#$ -cwd
##+++++++++++++++++++++++++++++++++++++++++++++

#source ~/miniconda3/bin/activate
conda activate sra_tools

for SAMPLE_ID in `cat sra.txt`; do

cd ${SAMPLE_ID}

fasterq-dump --split-files ${SAMPLE_ID}.sra

cd ../

done 

#fasterq-dump will extract all the fastq files from the downloaded accessions.
#after this script runs the raw downloaded files can be deleted to save space.
