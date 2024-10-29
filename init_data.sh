#!/bin/bash

# Create the necessary directories
mkdir -p ref_seqs/ecoli
mkdir -p ref_seqs/bacteria
mkdir -p ref_seqs/fruit_fly

# Define the reference sequences and their target directories
declare -A ref_sequences
ref_sequences["NC_000913"]="ref_seqs/ecoli"
ref_sequences["NC_001133-48"]="ref_seqs/bacteria"
ref_sequences["NT_033777-79"]="ref_seqs/fruit_fly"
ref_sequences["NC_001224"]="ref_seqs/fruit_fly"
ref_sequences["NC_004353-54"]="ref_seqs/fruit_fly"
ref_sequences["NT_037436"]="ref_seqs/fruit_fly"

# Download each reference sequence
for ref_id in "${!ref_sequences[@]}"; do
    target_dir="${ref_sequences[$ref_id]}"
    echo "Downloading FASTA file for reference sequence ID: $ref_id into $target_dir"
    esearch -db nucleotide -query "$ref_id" | efetch -format fasta > "$target_dir/${ref_id}.fasta"
done

echo "All reference sequences have been downloaded to their respective folders."
