# Chromatin Workshop  12/8/2024

## CSHL Fall 2024 Computational Genomics Course
Welcome! This is the material and tutorials for the Chromatin Workshop.
The database adopted in this course is under the reference: "Enhancer Chromatin and 3D Genome Architecture Changes from Naive to Primed Human Embryonic Stem Cell States". https://www.sciencedirect.com/science/article/pii/S2213671119301262?via%3Dihub <br />

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
- deepTools: https://deeptools.readthedocs.io/en/develop  <br />
- Samtools: https://www.htslib.org/  <br />
- Chromap: https://github.com/haowenz/chromap <br />
- bedtools: https://bedtools.readthedocs.io/en/latest/  <br />
- MACS2: https://pypi.org/project/MACS2/  <br />

## 0) Interactive Session and file setup
- The data is miniaturized, so we can run it in real time using an interactive session.
- Scripts are still provided for all of these in the directory

- This works the same way as submitting a script, but in one line and with qlogin (for AGE)
`qlogin -l mfree=5G -pe threads 8`
- Now we will copy the files from the shared directory to your own directory to play with.
`cp -r ~/../shared/CSHL_Chromatin_Workshop_2024 ~`


## Install deepTools & MultiQC <br />

***Do not run!!***
/grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
**0.1)** Install miniconda in your home directory <br />
0.1.1) `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh` *remember to change the permissions* <br />
&#x1F538; check here for others versions: https://docs.conda.io/projects/miniconda/en/latest/ <br />
 0.1.2) `bash Miniconda3-latest-Linux-x86_64.sh` *install miniconda* <br />
 0.1.3) activate base: `source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate` *check if is installed conda list* <br />
 **0.2)** Install deepTools <br />
 0.2.1) `conda create --name deeptools` *create a new environment* <br />
 0.2.2) `conda activate deeptools` *activate new environment* <br />
 0.2.3) `conda install -c bioconda deeptools` *download packages* <br /> 
 **0.3)** Install MultiQC <br /> 
 0.3.1) `conda create --name multiqc` *create a new environment* <br />
 0.3.2) `conda activate fastqc multiqc` *activate new environment* <br />
 0.3.3) `conda install -c bioconda multiqc` *download packages* <br /> 
 **0.4)** Install chromap <br /> 
 0.4.1) `conda create --name chromap` *create a new environment* <br />
 0.4.2) `conda activate chromap` *activate new environment* <br />
 0.4.3) `conda install -c bioconda chromap` *download packages* <br /> 
 **0.5)** Install macs <br /> 
 0.5.1) `conda create --name macs2` *create a new environment* <br />
 0.5.2) `conda activate macs2` *activate new environment* <br />
 0.5.3) `conda install -c bioconda macs2` *download packages* <br /> 
 **0.6)** Install basic_tools <br /> 
 0.6.1) `conda create --name basic_tools` *create a new environment* <br />
 0.6.2) `conda activate basic_tools` *activate new environment* <br />
 0.6.3) `conda install -c bioconda bedtools samtools seqkit dos2unix` *download packages* <br /> 

**0.7)** *OPTIONAL* Install sra-tools and seqtk <br />
- *you only need this if you want to prep the data yourself, but don't do it for this workshop*
 0.7.1) `conda create --name sra_tools` *create a new environment* <br />
 0.7.2) `conda activate sra_tools` *activate new environment* <br />
 0.7.3) `conda install -c bioconda sra-tools seqtk` *download packages* <br /> 



## 1) Quality Control of Sequencing using FastQC/MultiQC
- **Documentation**: *FastQC*: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ & *MultiQC*: https://multiqc.info/ <br />
`conda activate multiqc` <br />
`cd /grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024/data/subset/` <br />
- FastQC: ~3 -4 min <br /> 
`fastqc *.fastq`
- MultiQC: ~ 1 min <br />
`multiqc *.zip`

### 2) Processing data & Genome Mapping
**Chromap** https://github.com/haowenz/chromap for aligning and preprocessing high throughput chromatin profiles (*ATAC-seq & ChIP-seq*):
- **1)** Trim the low-quality reads and adapters
- **2)** Remove duplicated reads
- **3)** Perform the mapping. <br />

- Build the indexed genome ~ 1 min
```bash
conda activate chromap
cd /grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024/genome/
chromap -i -r hg38_chr22.fa -o index
```
- Flags:
**-i**: indexing genome flag 
**-r**: reference genome
**-o**: output name

### 2.1) Reference Genome
- Genome assembly: GRCh38.p14 <br />
&#x1F538; **Only the chromosome 22 human genome** (*only because it is one of the smallest!!*)
- Chromosome 22 info : https://useast.ensembl.org/Homo_sapiens/Location/Chromosome?r=22 <br />

### 2.2) Processing & mapping data 
- Chromap: *~35 sec* <br />
Chromap performs the remove duplicates, adapters and alignment using high throughput chromatin profiles. <br />
***If you use a different aligner you will need to do those steps yourself.***
```bash
#THESE ARE NOT ABSOLUTE PATHS YOU NEED TO ADD TO THEM
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir="${datadir}/data/subset"
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt

#source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate chromap

cd $dir
for SAMPLE_ID in `cat $list`; do
# map using chromap, output is bed file
chromap --preset chip -x $index -r $genome -q 20 --min-read-length 10   -1  ${SAMPLE_ID}.fastq  -o  ${SAMPLE_ID}_chromap.bed
done
conda activate basic_tools
dos2unix *.bed
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
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
cd ${datadir}/genome
conda activate basic_tools
samtools faidx hg38_chr22.fa
cut -f1,2 hg38_chr22.fa.fai > sizes.genome
cat sizes.genome
``` 

Now you can actually convert to the bam files for peak calling.
```bash
#THESE ARE INCOMPLETE PATHS, YOU NEED TO ADD TO THEIR BEGINNING
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir=${datadir}/data/subset
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt
size=${datadir}/genome/sizes.genome #size of chrm
##

cd $dir

# source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate basic_tools

for SAMPLE_ID in `cat $list`; do
#convert alignments to BAM
        bedtools bedtobam  -i ${SAMPLE_ID}_chromap.bed -g $size > ${SAMPLE_ID}_chromap.bam #input files are the same from chromap
done
```


2.3.2) Sort .**bam** & index generation **.bai** & convert to ***.bw** (*Big*Wig*) *~3min*  <br />
```bash
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
  samtools sort ${SAMPLE_ID}_chromap.bam -@ ${NSLOTS}  -o ${SAMPLE_ID}_chromap_sorted.bam
##remove
#samtools index ${SAMPLE_ID}_CT_chromap_sorted.bam
#samtools idxstats ${SAMPLE_ID}_CT_chromap_sorted.bam | cut -f1 | grep -v Mt | xargs samtools view -b ${SAMPLE_ID}_chromap_sorted.bam > ${SAMPLE_ID$
##sort
  samtools sort ${SAMPLE_ID}_chromap_sorted.bam -@ ${NSLOTS}  -o ${SAMPLE_ID}_treat.bam
##convert to bw
  samtools index ${SAMPLE_ID}_treat.bam
  conda activate deeptools
  bamCoverage -p max -b ${SAMPLE_ID}_treat.bam  --normalizeUsing RPKM  -v  -o ${SAMPLE_ID}_norm.bw ##you can use this on the genome browser
  conda activate basic_tools
done
```
create an index, convert bam to bw, & normalize data RPKM (deeptools) <br />

- **Extra** Remove the Chrm MT; &#x1F538; Should read the Mt in the reference genome *(check the reference and annotation genome)*; idxstats: create the index. <br />
Chromosome MT (Mitochondrial) can cause noise in the *calling peaks* should remove from the *.bam files  <br />
*DO NOT RUN* <br />
```
## Hi don't run me.
samtools index ${sorted.bam.file} 
samtools idxstats ${sorted.bam.file} | cut -f1 | grep -v Mt | xargs samtools view -b ${sorted.bam.file}  > ${sorted-noMT.bam.file}
```

 **check files**: Output file (*BAM format*); 
 - *Check the bigwig files in the genome browser*  <br />
`samtools view SRR5063143_naive_H3K27ac_chromap.bam | head -n 5` <br />
- Use the IGV app: https://igv.org/app/  <br />
- Select the hg38 genome and select the chromosome 22 *(chr22)* <br />
- download the *.bw files from the HPC to your personal PC *can use SCP or STFP*  <br />
***Whatever method to get files onto your computer is fine*** <br />
- upload the *.bw in the *track* function  <br />
- Let's have fun!! check the **FBXO7** gene  <br />
  *https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000100225;r=22:32474676-32498829
  *chr22:32474676-32498829   <br />
 
## 3) Peak Calling 
**MACS2** the Model-based Analysis of ChIP-Seq (MACS) for chormatin data analysis https://pypi.org/project/MACS2/ <br />
**Analysis for ChIP-seq; ATAC-seq; Cut&Tag**. *The parameters depend on the data type.*  <br />

- Histone modification" *H3K27ac:* broad peaks; *H3K4me3* narrow peaks. <br />   

macs2 *~2 min* <br />  

```bash
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

IFS=',' read -r file1 file2 <<< $files
##macs2-broad
macs2 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g hs --nomodel --shift -100 --extsize 200 \
        -n ${file1%_treat.bam} --broad \
        --outdir . 2> ${file1%_treat.bam}_broad_macs2.log
done


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


```


##QC Analysis *~ 3 min* 
**Fraction of reads in peaks (FRiP):** FRiP Score essential to evaluate the Peaks Quality. *more details:* https://yiweiniu.github.io/blog/2019/03/Calculate-FRiP-score/ <br />
- Request data: *.bam files & *Peaks files (Narrow or broad)

```bash

conda activate basic_tools
for peakfile in `ls *.Peak`; do

readfile=${peakfile//_peaks.Peak/_chromap_sorted.bam}

reads=$(samtools view -c ${readfile})
reads_peaks=$(bedtools intersect -u -a ${readfile} -b ${peakfile} -ubam | samtools view -c)

frip_score=$(echo "scale = 6; ${reads_peaks} / ${reads}" | bc)

echo -e "${peakfile//_peaks.Peak/}\t${frip_score}"
done

``` 

These FRiP values are terrible! Usually we would want a FRiP score of at least .2 or so. <br />
However, these are very subsampled, and only one chromosome is present, so it will be good for now.


## Let's visualize our results
**1)** deeptools:Correlation and Heatmap plots. Correlation matrix bewteen the replicates (*QC analysis*) and Heatmap (*visualize the signal intensity:Input; HK3me4;HK27ac*)  *~ 20 min* <br /> 

```bash
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir=${datadir}/data/subset

cd ${dir}
mkdir -p consensus_matrixes/
mkdir -p deeptools_graphs/

#source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate deeptools

multiBigwigSummary bins -b *norm.bw -o bw_corr.npz -p ${NSLOTS}

plotCorrelation -in bw_corr.npz -c spearman -p heatmap -o deeptools_graphs/correlation_heatmap.pdf
#plotCorrelation -in bw_corr.npz -c spearman -p scatterplot -o deeptools_graphs/correlation_scatterplot.pdf

```
Now to visualize the peak files

```bash
datadir="/grid/genomicscourse/home/beuchat/CSHL_Chromatin_Workshop_2024"
dir=${datadir}/data/subset

cd ${dir}

mkdir -p consensus_matrixes/
mkdir -p deeptools_graphs/

#source /grid/genomicscourse/home/shared/conda/miniconda3/bin/activate
conda activate deeptools

input_combinations=(
"SRR5063143_naive_H3K27ac_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063144_naive_H3K27ac_treat.bam,SRR5063154_naive_input_treat.bam"
"SRR5063149_naive_H3K4me3_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063150_naive_H3K4me3_treat.bam,SRR5063154_naive_input_treat.bam"
)

for files in ${input_combinations[@]}; do

IFS=',' read -r file1 file2 <<< $files

bamCompare -b1 ${file1} -b2 ${file2} -o ${file1//_treat\.bam/_differential}.bw
peakfile=${file1//_treat.bam/_peaks\.Peak}
bw_file=$(echo ${file1//_treat\.bam/_differential}.bw)

computeMatrix scale-regions -p ${NSLOTS} -S ${bw_file} -R ${peakfile} -b 3000 -a 3000 \
        -o consensus_matrixes/${peakfile//_peaks\.Peak/\.matrix}

plotHeatmap -m consensus_matrixes/${peakfile//_peaks\.Peak/\.matrix} -o deeptools_graphs/${peakfile//_peaks\.Peak/\.pdf} \
        --dpi 300 --startLabel "Peak Start" --endLabel "Peak End" -x "Distance" --heatmapWidth 12 --regionsLabel "Peaks"

done

```

## Bedtools! <br />

Bedtools is the classic way to do presence/absence analysis. <br />
More quantitative methods are available, such as through diffbind. <br />
However, presence/absence is still great to get a good idea about your data. <br />

```bash
conda activate basic_tools
#how many peaks do we have?
wc -l SRR5063143_naive_H3K27ac_peaks.Peak
wc -l SRR5063149_naive_H3K4me3_peaks.Peak

#the following command keeps only peaks present in both files.
bedtools intersect -a SRR5063143_naive_H3K27ac_peaks.Peak -b SRR5063149_naive_H3K4me3_peaks.Peak > out.bed
wc -l out.bed
#the following command keeps only peaks present in file 1 but not 2
bedtools intersect -v -a SRR5063143_naive_H3K27ac_peaks.Peak -b SRR5063149_naive_H3K4me3_peaks.Peak > out.bed
wc -l out.bed
#the following command keeps only peaks present in file 2 but not 1
bedtools intersect -b SRR5063143_naive_H3K27ac_peaks.Peak -a SRR5063149_naive_H3K4me3_peaks.Peak > out.bed
wc -l out.bed

```

##

