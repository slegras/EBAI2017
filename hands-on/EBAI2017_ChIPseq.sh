#! /bin/env bash

###################################################
################# Variable definition
login=slegras
home=/shared/projects/training/${login}/EBA2017_chipseq

###################################################
################# Working environment
cd /shared/projects/training/${login}

## Create a working directory for entire hands-on part
mkdir EBA2017_chipseq
cd EBA2017_chipseq

## Create a directory for raw data
mkdir data
cd data
## copy fastq/fasta files from local computer
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

## Run fastqc on IP
srun fastqc ../data/SRR576933.fastq.gz -o .

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
srun bowtie ../index/Escherichia_coli_K12 ../../data/SRR576933.fastq -v 2 -m 1 -3 1 -S 2> SRR576933.out > SRR576933.sam

## Compress back fastq IP file
srun gzip ../../data/SRR576933.fastq

## Go back to upper directory
cd ..

## Create a directory for IP alignment
mkdir Control

## Go to Control directory
cd Control

## Unzip fastq IP file
srun gunzip ../../data/SRR576938.fastq.gz

## Run alignment
srun bowtie ../index/Escherichia_coli_K12 ../../data/SRR576938.fastq -v 2 -m 1 -3 1 -S 2> SRR576938.out > SRR576938.sam

## Compress back fastq IP file
srun gzip ../../data/SRR576938.fastq

## Go to home working directory
cd ../..
