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

multiBigwigSummary bins -b *norm.bw -o bw_corr.npz -p ${NSLOTS}
#the above command takes all the bigwigs listed after -b separated by a space (automatically done here
    #with the *norm.bw), splits them into bins, figures out the pairwise correlations for each bin for each bigwig
    #averages the correlation among all the bin pairwise comparisons for each bigwig (so you get one correlation
    #value for each bigwig to bigwig comparison) and stores that in the output matrix.
    #-p ${NSLOTS} is for the number of parallel threads available to speed this calculation up.

plotCorrelation -in bw_corr.npz -c spearman -p heatmap -o deeptools_graphs/correlation_heatmap.pdf
#plotCorrelation -in bw_corr.npz -c spearman -p scatterplot -o deeptools_graphs/correlation_scatterplot.pdf

#This just plots a heatmap and clusters them hierarchicaly.
#the commented out scatterplot also works and essentially plots the same thing, but I find it less useful.