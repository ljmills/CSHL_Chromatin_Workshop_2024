#!/bin/bash
#$ -S /bin/bash
#$ -pe threads 1
#$ -l m_mem_free=80G -l h_rt=170:00:00
#$ -cwd
##+++++++++++++++++++++++++++++++++++++++++++++

# source ~miniconda3/bin/activate
conda activate sra_tools
prefetch --option-file sra.txt

#this downloads the list of accessions in the sra.txt document from SRA, which is where a lot of
#data that gets published in papers gets deposited.
#the prefetch command will download a compressed nonsnensical version of the file which must then
#be extracted with the fasterq-dump function. To get at that function, check the split.sh script.
