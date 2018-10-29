#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -m abe
#PBS -j oe
#PBS -M jholmes5@illinois.edu
#PBS -N trim_arr
#PBS -t 1-4

module load trimmomatic/0.33

cd $PBS_O_WORKDIR

line=$(sed -n -e "$PBS_ARRAYID p" brooke_basenames.txt)

java -jar /home/apps/trimmomatic/trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 4 -phred33  \
../data/raw-seq/${line}R1_001.fastq ../data/raw-seq/${line}R2_001.fastq ../results/trimmomatic/${line}R1_trim_paired.fastq \
../results/trimmomatic/${line}R1_trim_unpaired.fastq ../results/trimmomatic/${line}R2_trim_paired.fastq \
../results/trimmomatic/${line}R2_trim_unpaired.fastq ILLUMINACLIP:/home/apps/trimmomatic/trimmomatic-0.33/adapters/TruSeq3-PE-2.fa:2:15:10 MINLEN:30

