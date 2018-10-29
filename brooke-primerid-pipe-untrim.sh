#!/bin/bash
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=100GB
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -J primerid-pipe
#SBATCH -a 1-4

cd /home/groups/hpcbio/projects/brooke/2016-10-software/ 
line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" src/brooke_basenames_short.txt)


## Step 3 ##
 # Align with bwa mem and get coverage to see how many amplicons to expect 
 # and whether a gap or not.

echo "Starting Step 3"
 
# I did not update to the latest versions in case it breaks something, but you can try updating these 
# and see if it still works. 
module load BWA/0.7.17-IGB-gcc-4.9.4
module load SAMtools/1.4.1-IGB-gcc-4.9.4 
module load Java/1.8.0_121  # Previously this was run with version 65, but this is the oldest available on new Biocluster.
module load R/3.2.5-IGB-gcc-4.9.4

# makes index if it doesn't exist
if [ ! -f data/genome/HA_PR8.fa.bwt ]; then
	bwa index data/genome/HA_PR8.fa
fi

bwa mem -t $SLURM_CPUS_PER_TASK data/genome/HA_PR8.fa data/raw-seq/${line}_L001_R1_001.fastq \
data/raw-seq/${line}_L001_R2_001.fastq | samtools view -bS - > results/untrim-primerid/step3/${line}.bam \

java -Xmx2G -jar src/primer-id-progs/picard.jar SortSam I=results/untrim-primerid/step3/${line}.bam \
CREATE_INDEX=true O=results/untrim-primerid/step3/${line}.sort.bam SO=coordinate TMP_DIR=tmp/ \

graph_coverage.pl --bam_file results/untrim-primerid/step3/${line}.sort.bam --output_dir results/untrim-primerid/step3/ 

echo "Ending Step 3"

## Step 4 ##
 # Make contigs from overlapping reads.

echo "Starting Step 4"

module load PANDAseq/2.11-IGB-gcc-4.9.4

pandaseq -f data/raw-seq/${line}_L001_R1_001.fastq -r data/raw-seq/${line}_L001_R2_001.fastq \
-F -T $SLURM_CPUS_PER_TASK -u results/untrim-primerid/step4/${line}_pandaseq_unaligned.fastq > \
results/untrim-primerid/step4/${line}.contigs.fastq

module unload PANDAseq/2.11-IGB-gcc-4.9.4

echo "Ending Step 4"


## Step 5 ##
 # Extract primerID from sequences and place in the first line of the fastq record 

echo "Starting Step 5"

filter_fastq_by_primerid_length.pl --removepost --post CA --file_in \
results/untrim-primerid/step4/${line}.contigs.fastq --n 12 --output_dir \
results/untrim-primerid/step5/

echo "Ending Step 5"


## Step 6 ##
 # Split into amplicon regions and strip off primer sequences

echo "Starting Step 6"

Btrim64 -p src/2016-brooke-primers.txt -t results/untrim-primerid/step5/${line}.contigs.pid.fastq \
-o results/untrim-primerid/step6/${line}.contigs.pid.btrim.fastq -u 2 -v 2 -S -B -e 300

echo "Ending Step 6"


## Step 7 ##
 # Remove off-target sequences

echo "Starting Step 7"

module load BioPerl/1.7.1-IGB-gcc-4.9.4-Perl-5.24.1 

for fastq in results/untrim-primerid/step6/${line}.contigs.pid.btrim.fastq.*; \
do get_majority_block_bam.pl -ref data/genome/HA_PR8.fa -fastq $fastq \
 -output_dir results/untrim-primerid/step7/; done

echo "Ending Step 7"


## Step 8 ##
 # Optional step. Makes graph coverage (after cleaning up).

echo "Starting Step 8"

for bam in results/untrim-primerid/step7/${line}.contigs.pid.btrim.fastq.*bam; \
do graph_coverage.pl --bam_file $bam --output_dir results/untrim-primerid/step8/ ; done

echo "Ending Step 8"


## Step 9 ##
 # Merge the reads by primer id. First command creates *majority.group.counts.txt files.
 # Second command merges the reads.

#echo "Starting Step 9"

module load MAFFT/7.310-IGB-gcc-4.9.4

#Saves results inside of Step7/
for sample in results/untrim-primerid/step7/${line}.contigs.pid.btrim.*.majority.bam; \
do a=${sample/.majority.bam/}; merge_primerid_read_groups.pl \
--save results/untrim-primerid/step7/ -p 8 --plot_only --plot_counts $sample; \
done

echo "First command finished...Second command starting"

#Saves results inside of Step9/
for sample in results/untrim-primerid/step7/${line}.contigs.pid.btrim.*.majority.bam; 
do a=${sample/.majority.bam/}; groups=${a}.majority.group.counts.txt; 
cutoff=$(compute_cutoff.pl $(cat $groups | tail -n 1 | cut -f 1)); 
merge_primerid_read_groups.pl -m $cutoff --save results/untrim-primerid/step9/ --ambig 600 --min_freq 0.75 -p 8 $sample; 
done

echo "Ending Step 9"


## Step 10 ##
 # Converts merged reads to codons & amino acids to get frequency tables and clean 
 # read/peptide alignment files. This step can take a long time if you have more 
 # than 20,000 unique reads. See GitHub for alterantives.

echo "Starting Step 10"

cd results/untrim-primerid/step9/
for file in ${line}.contigs.pid.btrim.*.majority.cons.fasta; \
do convert_reads_to_amino_acid.pl -p 10 --files $file --ref ../../../data/genome/HA_PR8.fa --prefix ../step10/${file/.fasta/}; \
done 

echo "Ending Step 10"


## Step 11 ##
 # Calculates LD for all variants above a certain frequency level

echo "Starting Step 11"

cd ../step10/
for v in 0.005; do for i in ${line}; do for n in 0 1 2 3; \
do sample=${i}.contigs.pid.btrim.fastq.${n}.majority.cons.variants.minfreq0.xls; \
calculate_linkage_disequilibrium.pl $sample ${sample/.variants.minfreq0.xls/}.cleanreads.txt \
${sample/.variants.minfreq0.xls/}.cleanpeptides.txt --variant_threshold $v  \
--group_id Amplicon-${n} --label Sample-${i} --save ../step11/; done; done; done

echo "Ending Step 11"





