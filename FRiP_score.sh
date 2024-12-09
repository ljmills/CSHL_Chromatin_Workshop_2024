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


conda activate basic_tools
for peakfile in `ls *.Peak`; do
#the command in the tickmarks this time lists every file that ends with .Peak

readfile=${peakfile//_peaks.Peak/_chromap_sorted.bam}
#this line calls the variable $peakfile you just made, finds _peaks.Peak within it,
    #and replaces it with _chromap_sorted.bam,
    #which are the reads that made up the peak file.

reads=$(samtools view -c ${readfile})
#samtools view -c just counts the number of reads in the sam/bam file
reads_peaks=$(bedtools intersect -u -a ${readfile} -b ${peakfile} -ubam | samtools view -c)
#the bedtools intersect tells you the reads in -a that intersect with the peaks in -b,
    #the -u and -ubam flag means it should be printed as an unsorted bam,
    #then the samtools view -c counts THOSE to get the reads in peaks
frip_score=$(echo "scale = 6; ${reads_peaks} / ${reads}" | bc)
#to get the frip score we just divide one by the other. The only complication is that
    #bash does not deal with floating point numbers (numbers that have decimals on them)
    #so a simple $reads / $reads_peaks won't work. Instead we need to pass the expression through
    #a pipe to bc, which can. we do this by building our expression and printing it with echo.
    #scale = 6 is saying we should calculate and print to 6 decimal places.

echo -e "${peakfile//_peaks.Peak/}\t${frip_score}"
#this line just prints the name of the library (so peakfile but without the _peaks.Peak ending)
    #and the frip score separated by a tab.
done
