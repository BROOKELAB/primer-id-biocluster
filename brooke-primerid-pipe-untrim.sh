#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l mem=100GB
#PBS -m abe
#PBS -j oe
#PBS -M jholmes5@illinois.edu
#PBS -N primerid-pipe
#PBS -t 1-4

cd /home/groups/hpcbio/projects/brooke/2016-10-software/ 
line=$(sed -n -e "$PBS_ARRAYID p" src/brooke_basenames_short.txt)


## Step 3 ##
 # Align with bwa mem and get coverage to see how many amplicons to expect 
 # and whether a gap or not.

#echo "Starting Step 3"
 
module load bwa/0.7.15
module load samtools/1.3.1
module load java/1.8.0_65
module load R/3.2.3

# makes index if it doesn't exist
if [ ! -f data/genome/HA_PR8.fa.bwt ]; then
	bwa index data/genome/HA_PR8.fa
fi

#bwa mem -t $PBS_NUM_PPN data/genome/HA_PR8.fa data/raw-seq/${line}_L001_R1_001.fastq \
#data/raw-seq/${line}_L001_R2_001.fastq | samtools view -bS - > results/untrim-primerid/step3/${line}.bam \

#java -Xmx2G -jar src/primer-id-progs/picard.jar SortSam I=results/untrim-primerid/step3/${line}.bam \
#CREATE_INDEX=true O=results/untrim-primerid/step3/${line}.sort.bam SO=coordinate TMP_DIR=tmp/ \

#graph_coverage.pl --bam_file results/untrim-primerid/step3/${line}.sort.bam --output_dir results/untrim-primerid/step3/ 

#echo "Ending Step 3"

## Step 4 ##
 # Make contigs from overlapping reads.

#echo "Starting Step 4"

module load pandaseq/2.10

#pandaseq -f data/raw-seq/${line}_L001_R1_001.fastq -r data/raw-seq/${line}_L001_R2_001.fastq \
#-F -T $PBS_NUM_PPN -u results/untrim-primerid/step4/${line}_pandaseq_unaligned.fastq > \
#results/untrim-primerid/step4/${line}.contigs.fastq

module unload pandaseq/2.10

#echo "Ending Step 4"


## Step 5 ##
 # Extract primerID from sequences and place in the first line of the fastq record 

#echo "Starting Step 5"

#filter_fastq_by_primerid_length.pl --removepost --post CA --file_in \
#results/untrim-primerid/step4/${line}.contigs.fastq --n 12 --output_dir \
#results/untrim-primerid/step5/

#echo "Ending Step 5"


## Step 6 ##
 # Split into amplicon regions and strip off primer sequences

#echo "Starting Step 6"

#Btrim64 -p src/2016-brooke-primers.txt -t results/untrim-primerid/step5/${line}.contigs.pid.fastq \
#-o results/untrim-primerid/step6/${line}.contigs.pid.btrim.fastq -u 2 -v 2 -S -B -e 300

#echo "Ending Step 6"


## Step 7 ##
 # Remove off-target sequences## Step 7 ##

#echo "Starting Step 7"

module load perl/5.16.1_hpcbio
module load bioperl/bioperl-live 

#for fastq in results/untrim-primerid/step6/${line}.contigs.pid.btrim.fastq.*; \
#do get_majority_block_bam.pl -ref data/genome/HA_PR8.fa -fastq $fastq \
# -output_dir results/untrim-primerid/step7/; done

#echo "Ending Step 7"


## Step 8 ##
 # Optional step. Makes graph coverage (after cleaning up).

#echo "Starting Step 8"

#for bam in results/untrim-primerid/step7/${line}.contigs.pid.btrim.fastq.*bam; \
#do graph_coverage.pl --bam_file $bam --output_dir results/untrim-primerid/step8/ ; done

#echo "Ending Step 8"


## Step 9 ##
 # Merge the reads by primer id. First command creates *majority.group.counts.txt files.
 # Second command merges the reads.

#echo "Starting Step 9"

module load mafft/7.130

#Saves results inside of Step7/
#for sample in results/untrim-primerid/step7/${line}.contigs.pid.btrim.*.majority.bam; \
#do a=${sample/.majority.bam/}; merge_primerid_read_groups.pl \
#--save results/untrim-primerid/step7/ -p 8 --plot_only --plot_counts $sample; \
#done

#echo "First command finished...Second command starting"

#Saves results inside of Step9/
#for sample in results/untrim-primerid/step7/${line}.contigs.pid.btrim.*.majority.bam; 
#do a=${sample/.majority.bam/}; groups=${a}.majority.group.counts.txt; 
#cutoff=$(compute_cutoff.pl $(cat $groups | tail -n 1 | cut -f 1)); 
#merge_primerid_read_groups.pl -m $cutoff --save results/untrim-primerid/step9/ --ambig 600 --min_freq 0.75 -p 8 $sample; 
#done

#echo "Ending Step 9"


## Step 10 ##
 # Converts merged reads to codons & amino acids to get frequency tables and clean 
 # read/peptide alignment files. This step can take a long time if you have more 
 # than 20,000 unique reads. See GitHub for alterantives.

#echo "Starting Step 10"

cd results/untrim-primerid/step9/
#for file in ${line}.contigs.pid.btrim.*.majority.cons.fasta; \
#do convert_reads_to_amino_acid.pl -p 10 --files $file --ref ../../../data/genome/HA_PR8.fa --prefix ../step10/${file/.fasta/}; \
#done 

#echo "Ending Step 10"


## Step 11 ##
 # Calculates LD for all variants above a certain frequency level

#echo "Starting Step 11"

cd ../step10/
for v in 0.005; do for i in ${line}; do for n in 0 1 2 3; \
do sample=${i}.contigs.pid.btrim.fastq.${n}.majority.cons.variants.minfreq0.xls; \
calculate_linkage_disequilibrium.pl $sample ${sample/.variants.minfreq0.xls/}.cleanreads.txt \
${sample/.variants.minfreq0.xls/}.cleanpeptides.txt --variant_threshold $v  \
--group_id Amplicon-${n} --label Sample-${i} --save ../step11/; done; done; done

#echo "Ending Step 11"





