#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -j oe
#PBS -M jholmes5@illinois.edu
#PBS -N fastqc
#PBS -t 1-4

module load fastqc/0.11.5

cd $PBS_O_WORKDIR

line=$(sed -n -e "$PBS_ARRAYID p" brooke_basenames.txt)

#Raw reads
fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}R1_001.fastq
fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}R2_001.fastq

#Trimmed reads


