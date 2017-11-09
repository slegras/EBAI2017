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
# srun samtools sort SRR576933.sam | samtools view -Sb > SRR576933.bam

## create an index for the bam file
#srun samtools index SRR576933.bam

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

## Go to home working directory
cd $home

###################################################
################# Bonus: PhantomPeakQualTools
## creating output directory for alignment
mkdir 00-PhantomPeakQualTools

## Go to newly created directory
cd 00-PhantomPeakQualTools

## Create a TagAlign file from the bam file
## No need (see documentation)
# srun samtools view -F 0x0204 -o - SRR576933.bam \
# | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"}
# else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' \
#   | gzip -c > SRR576933_experiment.tagAlign.gz

## Run phantompeakqualtools
Rscript ../scripts/phantompeakqualtools/run_spp.R -c=../02-Mapping/SRR576933.bam  -savp -out=SRR576933_IP_phantompeaks
