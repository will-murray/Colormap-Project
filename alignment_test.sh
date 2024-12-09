#usage : 

data_dir=$1
ref=$2



cd "$data_dir"
blasr lr.fasta "$ref_genome" --header > og.sam


# # blasr lr_corr.fasta "$ref_genome" --header > ../blasr_alignments/og_.sam
