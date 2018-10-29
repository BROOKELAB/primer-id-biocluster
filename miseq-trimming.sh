#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -J trim_arr
#SBATCH -a 1-4

module load Trimmomatic/0.38-Java-1.8.0_152

cd $SLURM_SUBMIT_DIR

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" brooke_basenames.txt)

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 4 -phred33  \
../data/raw-seq/${line}R1_001.fastq ../data/raw-seq/${line}R2_001.fastq ../results/trimmomatic/${line}R1_trim_paired.fastq \
../results/trimmomatic/${line}R1_trim_unpaired.fastq ../results/trimmomatic/${line}R2_trim_paired.fastq \
../results/trimmomatic/${line}R2_trim_unpaired.fastq ILLUMINACLIP:/home/apps/trimmomatic/trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:15:10 MINLEN:30

