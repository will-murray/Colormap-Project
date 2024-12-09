#!/bin/bash


# Usage
# provide 4 arguments:
#   1. output dir
#   2. long read file name (.fasta)
#   3 & 4: paired-end short reads file names (.fastq)
#
# Example for test data:
# ./short_long_align.sh test_data pac.fasta ill_1.fastq ill_2.fastq
# ./short_long_align.sh ecoli_data/ SRR10971019_sub.fasta SRR13921543_1_sub.fastq SRR13921543_2_sub.fastq
#



# awk 'NR <= 50000 {print}' ecoli_data/SRR13921543_1.fastq > ecoli_data/SRR13921543_1_sub.fastq
# awk 'NR <= 50000 {print}' ecoli_data/SRR13921543_2.fastq > ecoli_data/SRR13921543_2_sub.fastq
# awk 'NR <= 10000 {print}' ecoli_data/SRR10971019.fasta > ecoli_data/SRR10971019_sub.fasta



# Exit on error
set -e

# Arguments
output_dir=$1
long_read_fasta=$2
short_read_fastq_1=$3
short_read_fastq_2=$4

# Change to the output directory
cd "$output_dir"
bwa index "$long_read_fasta"
bwa mem -t 4 -aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 10 -w 40 -r1 -D 0 -y 20 -L 30,30 -T 2.5 "$long_read_fasta" "$short_read_fastq_1" "$short_read_fastq_2" > sl_align.sam

# Extract relevant alignments
grep -v "^@" sl_align.sam | #remove filter for non header lines from the same file
awk '{print $1 "." NR%2, $3, $4, $4 + length($10), $10}' |
sort -k2,2 -k3,3 -k4,4 | awk '$5 != "*" {print}' > sl_raw_align.txt



