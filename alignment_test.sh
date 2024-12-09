#usage : 

data_dir=$1
ref=$2



cd ecoli
blasr lr.fasta NC_000913.fasta   --header > og.sam
blasr lr_corr.fasta NC_000913.fasta --header > corr.sam