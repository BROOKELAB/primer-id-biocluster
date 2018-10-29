#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -J fastqc
#SBATCH -a 1-4

module load FastQC/0.11.5-IGB-gcc-4.9.4-Java-1.8.0_152

cd $SLURM_SUBMIT_DIR

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" brooke_basenames.txt)

#Raw reads
fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}R1_001.fastq
fastqc -outdir=../results/fastqc_raw/ --noextract ../data/raw-seq/${line}R2_001.fastq
