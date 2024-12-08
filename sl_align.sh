#!/bin/bash

# Usage:
# Provide 4 arguments:
#   1. output dir
#   2. long read file name (.fasta)
#   3. paired-end short read file 1 (.fastq)
#   4. paired-end short read file 2 (.fastq)
# Example:
# ./short_long_align.sh test_data pac.fasta ill_1.fastq ill_2.fastq

# Exit on error
set -e

# Arguments
output_dir=$1
long_read_fasta=$2
short_read_fastq_1=$3
short_read_fastq_2=$4
chunk_size=50000  # Number of reads per chunk
k=$5  # Number of chunks to process for testing

# Change to the output directory
cd "$output_dir"
mkdir "fastq_chunks"