#usage : 

data_dir=$1
og_long_reads=$2      #origonal long reads
corr_long_reads=$3    #corrected longs reads
ref_genome=$4         #reference genome


cd "$data_dir"
blasr "$og_long_reads" "$ref_genome" --header > ../blasr_alignments/og.sam

# blasr $corr_long_reads $refgenome  -noSplitSubreads -bestn 1 --header > ../blasr_alignments/corr.sam
# blasr SRR10971019_corr.fasta NC_000913.fasta