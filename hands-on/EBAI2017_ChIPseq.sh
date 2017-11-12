#! /bin/env bash

###################################################
################# Variable definition
login=slegras
home=/shared/projects/training/${login}/EBA2017_chipseq

###################################################
################# Working environment
cd ${home}

## Create a working directory for entire hands-on part
mkdir EBA2017_chipseq
cd EBA2017_chipseq

## Create a directory for raw data
mkdir data
cd data
## copy fastq/fasta files from local computer
## fastq file got downloaded from EBI
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576933/SRR576933.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576938/SRR576938.fastq.gz
## genome fasta file got downloaded from NCBI website
## annotation track (.gtf) was downloaded from UCSC table browser (gtf file)
## annotation track (tsv) was downloaded from UCSC table browser (selecting fields to be output)
# http://microbes.ucsc.edu/cgi-bin/hgTables?org=Escherichia+coli+K12&db=eschColi_K12&hgsid=1465191&hgta_doMainPage=1

## Changing chromosome names
srun zcat Escherichia_coli_K_12_MG1655.annotation.gtf | perl -pe 's/^chr/gi\|49175990\|ref\|NC_000913.2\|/' | \
 gzip > Escherichia_coli_K_12_MG1655.annotation.fixed.gtf.gz

cd $home

## Create a directory for scripts
mkdir scripts

## Go to the newly created directory
cd scripts

## Clone phantompeakqualtools repository
git clone https://github.com/crazyhottommy/phantompeakqualtools.git

## go back to home working directory
cd $home

## Loading conda ChIP-Seq environment
source activate eba2017_chipseq

###################################################
################# Quality controls
## Create a directory for quality control
mkdir 01-QualityControl

## Go to quality control directory
cd 01-QualityControl

## Test fastqc tool
## To be run with srun
srun fastqc --help

## Run fastqc on the IP
srun fastqc ../data/SRR576933.fastq.gz -o .

## Run fastqc on the control
srun fastqc ../data/SRR576938.fastq.gz -o .

## Go to home working directory
cd $home

###################################################
################# Mapping
## creating output directory for alignment
mkdir 02-Mapping

## Go to newly created directory
cd 02-Mapping

## Create a directory for index files
mkdir index

## Go to index directory
cd index

## testing bowtie-build
srun bowtie-build

## Unzip genome fasta file
srun gunzip ../../data/Escherichia_coli_K12.fasta.gz

## Creating genome index
srun bowtie-build ../../data/Escherichia_coli_K12.fasta Escherichia_coli_K12

## Compress back genome fasta file
srun gzip ../../data/Escherichia_coli_K12.fasta

## Go back to upper directory i.e 02-Mapping
cd ..

## Create a directory for IP alignment
mkdir IP

## Go to IP directory
cd IP

## Unzip fastq IP file
srun gunzip ../../data/SRR576933.fastq.gz

## Run alignment
sbatch bowtie ../index/Escherichia_coli_K12 ../../data/SRR576933.fastq -v 2 -m 1 -3 1 -S 2> SRR576933.out > SRR576933.sam

## Compress back fastq IP file
srun gzip ../../data/SRR576933.fastq

## Create a sorted bam file
srun samtools sort SRR576933.sam | samtools view -Sb > SRR576933.bam

## create an index for the bam file
srun samtools index SRR576933.bam

## Compress the .sam file
gzip SRR576933.sam

## Go back to upper directory
cd ..

## Create a directory for IP alignment
mkdir Control

## Go to Control directory
cd Control

## Unzip fastq IP file
srun gunzip ../../data/SRR576938.fastq.gz

## Run alignment
sbatch bowtie ../index/Escherichia_coli_K12 ../../data/SRR576938.fastq -v 2 -m 1 -3 1 -S 2> SRR576938.out > SRR576938.sam

## Compress back fastq IP file
srun gzip ../../data/SRR576938.fastq

## Create a sorted bam file
srun samtools sort SRR576938.sam | samtools view -Sb > SRR576938.bam

## create an index for the bam file
srun samtools index SRR576938.bam

## Compress the .sam file
gzip SRR576938.sam

## Go to the IP directory
cd ../IP

## Run picard markDuplicates on the IP sample
srun picard MarkDuplicates \
CREATE_INDEX=true \
INPUT=SRR576933.bam \
OUTPUT=Marked_SRR576933.bam \
METRICS_FILE=metrics \
VALIDATION_STRINGENCY=STRICT

## Go to the Control directory
cd ../Control

## Run picard markDuplicates on the Control sample
srun picard MarkDuplicates \
CREATE_INDEX=true \
INPUT=SRR576938.bam \
OUTPUT=Marked_SRR576938.bam \
METRICS_FILE=metrics \
VALIDATION_STRINGENCY=STRICT

## Go to home working directory
cd $home

###################################################
################# Bonus: PhantomPeakQualTools
## creating output directory for alignment
mkdir 00-PhantomPeakQualTools

## Go to newly created directory
cd 00-PhantomPeakQualTools

## create an R environment
# conda create -c r --name R r
# source activate R
# conda install -c bioconda r-spp
# conda install -c bioconda samtools
# conda install -c bioconda gawk

## convert the BAM file into TagAlign format, specific to the program that calculates the quality metrics
srun samtools view -F 0x0204 -o - ../02-Mapping/IP/SRR576933.bam | \
gawk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"}
else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' \
 | gzip -c > SRR576933_experiment.tagAlign.gz

## Run phantompeakqualtools
srun Rscript ../scripts/phantompeakqualtools/run_spp.R -c=SRR576933_experiment.tagAlign.gz  -savp -out=SRR576933_IP_phantompeaks

## Source back the chipseq environment
# source activate eba2017_chipseq

## Go to home working directory
cd $home

###################################################
################# visualization
## Create a directory to store visualization files
mkdir 03-Visualization

## Go to the newly created directory
cd 03-Visualization

## Test bamCoverage
srun bamCoverage --help

## Run bamCoverage on IP
srun --mem=3G bamCoverage --bam ../02-Mapping/IP/Marked_SRR576933.bam \
--outFileName SRR576933_nodup.bw --outFileFormat bigwig --normalizeTo1x 4639675 \
--skipNonCoveredRegions --extendReads 200 --ignoreDuplicates

## Run bamCoverage on Control
srun --mem=5G bamCoverage --bam ../02-Mapping/Control/Marked_SRR576938.bam \
--outFileName SRR576938_nodup.bw --outFileFormat bigwig --normalizeTo1x 4639675 \
--skipNonCoveredRegions --extendReads 200 --ignoreDuplicates

###################################################
################# Peak Calling
## Create a directory to store peak calling result files
mkdir 04-PeakCalling

## Go to the newly created directory
cd 04-PeakCalling

## Check macs parameters
srun macs

## Run macs on the IP and the Control file
srun macs -t ../02-Mapping/IP/SRR576933.bam -c ../02-Mapping/Control/SRR576938.bam --format BAM  --gsize 4639675 \
--name "FNR_Anaerobic_A" --bw 400 --bdg --single-profile --diag &> MACS.out

###################################################
################# Peak Annotation
## Create a directory to store peak annotation data
mkdir 05-PeakAnnotation

## Go to the newly created directory
cd 05-PeakAnnotation

## Uncompress annotation file
srun gunzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf.gz

## Create a BED file with 6 columns
awk -F "\t" '{print $0"\t+"}' ../04-PeakCalling/FNR_Anaerobic_A_peaks.bed > FNR_Anaerobic_A_peaks.bed

## Try it out
srun annotatePeaks.pl

## Annotation peaks with Homer
srun annotatePeaks.pl \
FNR_Anaerobic_A_peaks.bed \
../data/Escherichia_coli_K12.fasta \
-gtf ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf \
> FNR_Anaerobic_A_annotated_peaks.tsv

## Compress annotation file
srun gzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf

## Add gene symbol annotation using R
R
d <- read.table("FNR_Anaerobic_A_annotated_peaks.tsv", sep="\t", h=T)
gene.symbol <- read.table("../data/Escherichia_coli_K_12_MG1655.annotation.tsv.gz", h=F)
d.annot <- merge(d[,c(seq(1,6,1),8,10,11)], gene.symbol, by.x="Nearest.PromoterID", by.y="V1")
colnames(d.annot)[2] <- "PeakID"
colnames(d.annot)[dim(d.annot)[2]] <- "Gene.Symbol"
write.table(d.annot, "FNR_Anaerobic_A_final_peaks_annotation.tsv", col.names=T, row.names=F, sep="\t", quote=F)
quit()
n


###################################################
################# Motif analysis
## Create a directory to store data needed from motif analysis
mkdir 06-MotifAnalysis

## Go to the newly created directory
cd 06-MotifAnalysis

## Uncompress the genome file
srun gunzip ../data/Escherichia_coli_K12.fasta.gz

## Create an index of fasta file
srun samtools faidx ../data/Escherichia_coli_K12.fasta

## Extract fasta sequence from genomic coordinate of peaks
srun bedtools getfasta -fi ../data/Escherichia_coli_K12.fasta \
-bed ../04-PeakCalling/FNR_Anaerobic_A_peaks.bed -fo FNR_Anaerobic_A_peaks.fa

## Compress back the genome file
srun gzip ../data/Escherichia_coli_K12.fasta
