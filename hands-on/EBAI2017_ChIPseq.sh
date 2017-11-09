#! /bin/env bash

################# définition des variables
login=slegras
home=/shared/projects/training/${login}/EBA2017_chipseq

################# Gestion de l'environnement de travail
cd /shared/projects/training/${login}

## création d'un répertoire d'analyse
mkdir EBA2017_chipseq
cd EBA2017_chipseq

## création d'un répertoire qui va contenir les données brutes
mkdir data
cd data
## copie des fichiers fastq du TP depuis l'ordinateur local.
cd $home

## chargement de l'environnement conda pour le ChIP-Seq
source activate eba2017_chipseq

################# Lancement de l'analyse de la qualité des données
## creation du fichier de résultats des analyses de qualité
mkdir 01-QualityControl

## On va dans le répertoire que l'on vient de créer
cd 01-QualityControl

## test de lancement de fastqc
## On lance avec srun
srun fastqc --help

## Lancement de fastqc sur l'IP
srun fastqc ../data/SRR576933.fastq.gz -o .
