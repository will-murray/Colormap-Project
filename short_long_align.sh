# CLI's :
# 1. output dir
# 2. long read file name (.fasta)
# 3 & 4: paired end short reads file names (.fastq)

bwa index $1/$2
bwa mem $1/$2 $1/$4 $1/$5 > $1/sl_align.sam


