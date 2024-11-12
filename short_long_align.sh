#Usage
# provide 4 arguments:
#   1. output dir
#   2. long read file name (.fasta)
#   3 & 4: paired end short reads file names (.fastq)

# example for test data: ./short_long_align.sh test_data pac.fasta ill_1.fastq ill_2.fastq

cd $1
bwa index $2
-aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 16 -w 40 -r 1 -D 0 -y 20 -L 30,30 -T 2.5
bwa mem $2 $3 $4  > sl_align.sam


grep "^[ill]" sl_align.sam |
awk '{print $1 "." NR%2, $3, $4, $4 + length($10), $10}' |
sort -k2,2 -k3,3 -k4,4 > sl_raw_align.txt

awk '{count[$2]++} END {for (num in count) print num, count[num]}' sl_raw_align.txt > freq.txt





