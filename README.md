# Chromatin Workshop  02/21/2025

## CSHL Fall 2024 Computational Genomics Course
Welcome! This is the material and tutorials for the Chromatin Workshop.
The database adopted in this course is under the reference: "Enhancer Chromatin and 3D Genome Architecture Changes from Naive to Primed Human Embryonic Stem Cell States". https://www.sciencedirect.com/science/article/pii/S2213671119301262?via%3Dihub <br />
GEO page where the data is deposited: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90680 <br />

## Guidelines:
1) Setting up Conda <br />
2) QC Analysis & mapping data <br />
3) Pos-mapping data <br />
4) data visualization (IGV) <br />
5) Calling Peaks <br />
6) QC Peaks (FRiP Scores) <br />
7) Bedtools <br />

## Data info:
- **Cell Type**: Naïve cells <br />
- **Histone modification**: H3K4me3: *Promoters* & H3K27ac: *Promoters and Enhancers* <br />
- **Library info**: **1)** SE-fastq files; **2)** 3 M reads <br />
- **Data File**: Total 6 files (2 for each histone modification & 2 input files) <br />
##

## Tools and Packages Required: <br />
- miniconda (python)  https://docs.anaconda.com/miniconda/install/ <br />
- Sra-tools https://github.com/ncbi/sra-tools (module load sratoolkit/3.0.0)
- deepTools: https://deeptools.readthedocs.io/en/develop <br />
- Samtools: https://www.htslib.org/  <br />
- Chromap: https://github.com/haowenz/chromap <br />
- bedtools: https://bedtools.readthedocs.io/en/latest/ (module load bedtools2/2.31.0-gcc-8.2.0-7j35k74) <br />
- MACS2: https://pypi.org/project/MACS2/ (module load macs/2.1.1)  <br />
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ (module load fastqc/0.12.1)
- MultiQC: https://github.com/MultiQC/MultiQC

## 0) Interactive Session and Miniconda Setup 
- You will want to start an interactive session to do these installation. 
`srun --nodes=1 --ntasks-per-node=1 --mem=20g --time=4:00:00 -p agsmall --pty bash -i`
- Once that session is started we will install Miniconda . Download the script from the website (link above) into your MSI home directory (Linux x86) and run that script.
- create a chipseq conda enviroment
`conda create -n chipseq`
- activate this conda enviroment
  `conda activate chipseq`
  - You will be able to install the tools we need for this workshop inside of this conda enviroment. 


## 0.1)Install Software Not Available in Modules <br />
** Install deepTools <br />
 - `conda activate chipseq` *activate new environment if not already activated* <br />
 - `conda install -c conda-forge -c bioconda deeptools` *download packages* <br />
 ** Install Samtools
 - `conda install -c bioconda samtools` <br />
** Install MultiQC <br /> 
 - `conda install -c bioconda multiqc` *download packages* <br />
** Install chromap <br /> 
 - `conda install -c bioconda chromap` *download packages* <br /> 

** Install Other Tools <br /> 
 - `conda install -c bioconda seqkit dos2unix` *download packages* <br />
 - `conda install -c bioconda seqtk` *download packages* <br /> 

## 1) Download FASTQ files from GEO/SRA
- Create a directory in /scratch.global for this data `mkdir <yourWorkingDir>` then move into that directory `cd <yourWorkingDir>`
- Head to the GEO page, the SRA Run Selector at the bottom of the page.
- Select the FASTQ files you want to download and then create an Accession List. You will want to download the data for the H3K27ac, H3K27me3 and input from the Naive Cells. 6 samples total. 
- Put that accession list onto MSI (sftp, FileZilla, On Demand)
- `module load sratoolkit/3.0.0` this will give you access to sra-tools
- The first time you use sra-tools you will need to configure it
  `vdb-config -i` 
  `vdb-config --prefetch-to-cwd` This will tell prefetch to download files to the current working directory
- download fastq with prefetch then convert to fastq with fasterq-dump
- `cat SRR_Acc_List.txt | xargs prefetch` you should have 6 SRRXXXX directories
- `cat SRR_Acc_List.txt | xargs fasterq-dump` you should have 6 .fastq files
  
## 1) Quality Control of Sequencing using FastQC/MultiQC
- Run FastQC on all fastq files to look at the quality of the sequencing data.
- Make sure your chipseq conda environment is active and  you have loaded the fastqc module
- FastQC:  <br />
  `module load fastqc/0.12.1` 
  `fastqc *.fastq` 
- MultiQC: will combine all of your fastqc outputs into a single report <br />
`multiqc *.zip`

### 2) Processing data & Genome Mapping
**Chromap** https://github.com/haowenz/chromap for aligning and preprocessing high throughput chromatin profiles (*ATAC-seq & ChIP-seq*):
- **1)** Trim the low-quality reads and adapters
- **2)** Remove duplicated reads
- **3)** Perform the mapping. <br />

### 2.1) Reference Genome
- Genome assembly: GRCh38.p14 <br />
&#x1F538; **Only the chromosome 22 human genome** (*only because it is one of the smallest!!*)
- Chromosome 22 info : https://useast.ensembl.org/Homo_sapiens/Location/Chromosome?r=22 <br />
- You will need a reference genome to align you data too. And like all high throughput sequence alignment tools Chromap needs a specalized index for each reference genome you might want to align too.
- GRCh38 can be found here: `/common/bioref/ensembl/main/Homo_sapiens-113/GRCh38.p14`

- Build the indexed genome ~ 1 min. Only run the chr22 Only command
```bash
#chromap -i -r /common/bioref/ensembl/main/Homo_sapiens-113/GRCh38.p14/seq/genome.fa -o GRCh38_chromap_index
chromap -i -r /scratch.global/ljmills/GCD8141/chr22Only/Homo_sapiens.GRCh38.dna.chromosome.22.fa -o chr22Only_GRCh38_chromap_index

```
- Flags:
**-i**: indexing genome flag 
**-r**: reference genome
**-o**: output name



### 2.2) Processing & mapping data 
- Chromap: *~35 sec* <br />
Chromap performs the remove duplicates, adapters and alignment using high throughput chromatin profiles. <br />
***If you use a different aligner you will need to do those steps yourself.***
  Command for a single sample
```bash
chromap --preset chip -x chr22Only_GRCh38_chromap_index -r /scratch.global/ljmills/GCD8141/chr22Only/Homo_sapiens.GRCh38.dna.chromosome.22.fa -q 20 --min-read-length 10   -1  SRR5063143.fastq  -o  SRR5063143_chromap.bed

dos2unix SRR5063143_chromap.bed
```
- Flags:
**--preset chip**:Mapping chip reads
**-x**:index 
**-r**:reference genome 
**-q**:Min MAPQ (*Quality*)
**--min-read-length**:min length read
**-1**:Single-end (*include -2,implies they need both -1 and -2*) 
**-o**: output file
*--trim-adapters(not used)*

**check files**: Output file (*BED format*) <br /> 

`head SRR5063143_naive_H3K27ac_chromap.bed` <br />
- **chrm; start; end; N; q; strand**. <br />
22 &nbsp; 10510250 &nbsp; 10510300 &nbsp; N &nbsp; 59 &nbsp; + <br /> 
22 &nbsp; 10510252 &nbsp; 10510302 &nbsp; N &nbsp; 46 &nbsp; - <br /> 
22 &nbsp; 10511600 &nbsp; 10511650 &nbsp; N &nbsp; 60 &nbsp; + <br /> 

### 2.3) Post-mapping data 
2.3.1) Convert bed to bam *~2sec* <br /> 
- A note that these BAM files will be lacking many important components, but are usable for peak calling.

&#x1F538;**MUST!!** use the same version of reference genome use in the analysis <br />

Before you can convert to bams, you will need to calculate the size of each chromosome. <br />
We are only using chromosome 22, but the same commands will work with any genome.
```bash
samtools faidx --fai-idx chr22Only_GRCh38_faidx /scratch.global/ljmills/GCD8141/chr22Only/Homo_sapiens.GRCh38.dna.chromosome.22.fa 
cut -f1,2 chr22Only_GRCh38_faidx  > chr22sizes.genome
cat chr22sizes.genome
``` 

Now you can actually convert to the bam files for peak calling.
```bash
bedtools bedtobam  -i SRR5063143_chromap.bed -g chr22sizes.genome > SRR5063143_chromap.bam
```


2.3.2) Sort .**bam** & index generation **.bai** & convert to ***.bw** (*Big*Wig*) *~3min*  <br />
create an index, convert bam to bw, & normalize data RPKM (deeptools) <br />
```bash
samtools sort -o SRR5063143_chromap.sort.bam SRR5063143_chromap.bam
samtools index SRR5063143_chromap.sort.bam
###generate a bigWig file you can visulize in a genome browser 
bamCoverage -p max -b SRR5063143_chromap.sort.bam --normalizeUsing RPKM -v -o SRR5063143_chromap.norm.bw ##you can use this on the genome browser
```


## You ran all of this for 1 file but how do you effectively run all of these steps for all of your files?
You write a script and submit a job to MSI. I have a script started that you can edit to run these steps on all 6 of your FASTQ files. `cp /scratch.global/ljmills/GCD8141/classExampleScripts/chromatinWorkshop.sh <your directory>`
 
## 3) Peak Calling 
**MACS2** the Model-based Analysis of ChIP-Seq (MACS) for chormatin data analysis https://pypi.org/project/MACS2/ <br />
**Analysis for ChIP-seq; ATAC-seq; Cut&Tag**. *The parameters depend on the data type.*  <br />

- Histone modification" *H3K27me3:* broad peaks; *H3K4me3, H3K27ac* narrow peaks. <br />
- Which marks are broad vs narrow? https://www.encodeproject.org/chip-seq/histone/#histone
- You *NEED* a input control (-c) and your mark of interest (-t) to call peaks effectivly.
- Alignments for input controls can be found in `/scratch.global/ljmills/GCD8141/classAlignments`. copy them into your working directory

macs2 *~2 min* <br />  

```bash
# input 1
cp /scratch.global/ljmills/GCD8141/classAlignments/SRR5063153.chromap.sort.bam <your directory>
cp /scratch.global/ljmills/GCD8141/classAlignments/SRR5063154.chromap.sort.bam <your directory>

input1=SRR5063153.chromap.sort.bam
input2=SRR5063154.chromap.sort.bam

#H3K27ac narrow peak calling mode
macs2 callpeak \
-t SRR5063143.chromap.sort.bam \
-c $input1 \
-f BAM \
-g hs \
--nomodel --shift -100 \
--extsize 200 \
-n Naive_H3K27ac_rep1 \
--outdir <your directory>

#calculate FRiP (fraction of reads in peaks) 
readsInSample=$(samtools view -c SRR5063143.chromap.sort.bam)
readsInPeaks=$(bedtools intersect -u -a SRR5063143.chromap.sort.bam -b Naive_H3K27ac_rep1_peaks.narrowPeak -ubam | samtools view -c)
frip_score=$(echo "scale = 6; ${readsInPeaks} / ${readsInSample}" | bc)
echo -e "Naive_H3K27ac_rep1\t${frip_score}" >> FRiP_scores.tsv

#H3K27me3 broad peak calling mode
#H3K27me3
macs2 callpeak \
-t SRR5063145.chromap.sort.bam \
-c $input1 \
-f BAM \
-g hs \
--nomodel --shift -100 \
--extsize 200 \
-n Naive_H3K27me3_rep1 \
--broad \
--outdir <your dirctory> 

readsInSample=$(samtools view -c SRR5063145.chromap.sort.bam)
readsInPeaks=$(bedtools intersect -u -a SRR5063145.chromap.sort.bam -b Naive_H3K27me3_rep1_peaks.broadPeak -ubam | samtools view -c)
frip_score=$(echo "scale = 6; ${readsInPeaks} / ${readsInSample}" | bc)
echo -e "Naive_H3K27me3_rep1\t${frip_score}" >> FRiP_scores.tsv

```

- Again this is very laborous to write out for each input/treamtent pair. How do you keep track of what goes with what? While these peak calling jobs are fast it will take more time with the whole genome.
- Write a script to run these commands. Example script `/scratch.global/ljmills/GCD8141/classExampleScripts/callPeaks_FRiP.sh`

##QC Analysis *~ 3 min* 
**Fraction of reads in peaks (FRiP):** FRiP Score essential to evaluate the Peaks Quality. *more details:* https://yiweiniu.github.io/blog/2019/03/Calculate-FRiP-score/ <br />
- take a look at the FRiP_scores.tsv file you generated 

These FRiP values are terrible! Usually we would want a FRiP score of at least .2 or so. <br />
However, these are very subsampled, and only one chromosome is present, so it will be good for now.


## Let's visualize our results
**1)** deeptools:Correlation and Heatmap plots. Correlation matrix bewteen the replicates (*QC analysis*) and Heatmap (*visualize the signal intensity:Input; HK3me4;HK27ac*)  *~ 20 min* <br /> 

```bash

multiBigwigSummary bins -b *norm.bw -o bw_corr.npz -p ${NSLOTS}

plotCorrelation -in bw_corr.npz -c spearman -p heatmap -o correlation_heatmap.pdf
plotCorrelation -in bw_corr.npz -c spearman -p scatterplot -o correlation_scatterplot.pdf

```
Now to visualize the peak files

```bash
#H3K27ac
bamCompare -b1 SRR5063143.chromap.sort.bam \
-b2 SRR5063153.chromap.sort.bam \
-o H3K27ac_rep1_differential.bw

computeMatrix scale-regions  \
-S H3K27ac_rep1_differential.bw \
-R Naive_H3K27ac_rep1_peaks.narrowPeak \
 -b 3000 -a 3000 \
-o Naive_H3K27ac_rep1_peaks.matrix

plotHeatmap -m Naive_H3K27ac_rep1_peaks.matrix \
 -o Naive_H3K27ac_rep1_peaks.heatmap.pdf \
--dpi 300 \
--startLabel "Peak Start" \
--endLabel "Peak End" \
-x "Distance" \
--heatmapWidth 12 --regionsLabel "Peaks"

```

## Bedtools! <br />

Bedtools is the classic way to do presence/absence analysis. <br />
More quantitative methods are available, such as through diffbind. <br />
However, presence/absence is still great to get a good idea about your data. <br />

```bash
#how many peaks do we have?
wc -l Naive_H3K27ac_rep1_peaks.narrowPeak
wc -l Naive_H3K27me3_rep1_peaks.broadPeak

#the following command keeps only peaks present in both files.
bedtools intersect -a Naive_H3K27ac_rep1_peaks.narrowPeak \
-b Naive_H3K27me3_rep1_peaks.broadPeak > overlapingPeaks_Naive_H3K27ac_Naive_H3K27me3.bed

wc -l overlapingPeaks_Naive_H3K27ac_Naive_H3K27me3.bed

#the following command keeps only peaks present in file 1 but not 2
bedtools intersect -v -a Naive_H3K27ac_rep1_peaks.narrowPeak \
-b Naive_H3K27me3_rep1_peaks.broadPeak > Naive_H3K27ac_only.bed
wc -l Naive_H3K27ac_only.bed

#the following command keeps only peaks present in file 2 but not 1
bedtools intersect -v -a Naive_H3K27me3_rep1_peaks.broadPeak \
-b Naive_H3K27ac_rep1_peaks.narrowPeak > Naive_H3K27me3_only.bed
wc -l Naive_H3K27me3_only.bed

```

##

- **Extra** Remove the Chrm MT; &#x1F538; Should read the Mt in the reference genome *(check the reference and annotation genome)*; idxstats: create the index. <br />
Chromosome MT (Mitochondrial) can cause noise in the *calling peaks* should remove from the *.bam files  <br />
*DO NOT RUN* <br />
```
## Hi don't run me.
samtools index ${sorted.bam.file} 
samtools idxstats ${sorted.bam.file} | cut -f1 | grep -v Mt | xargs samtools view -b ${sorted.bam.file}  > ${sorted-noMT.bam.file}
```

