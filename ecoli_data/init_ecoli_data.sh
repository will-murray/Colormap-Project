# This script downloads the required files pertaining to Escherichia_coli_K12_MG1655, this includes:
# - reference genome : https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3/
# - PacBio long reads: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2063112&display=metadata
# - Illumina short reads: https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&page_size=10&acc=ERR022075&display=metadata


if [ ! -f SRR2063112.fasta ]; then
    echo "extracting SRR2063112 (long) with fastq-dump..."
    fastq-dump SRR2063112 
    echo "converting result to .fasta format..."
    awk '{if(NR%4==1) {print ">" substr($0, 2)} else if(NR%4==2) {print $0}}' SRR2063112.fastq > SRR2063112.fasta
    rm SRR2063112.fastq
fi


if [ ! -f ERR022075_1.fastq ]; then
    echo "extracting ERR022075 (short reads) with fasterq-dump..."
    fastq-dump --split-files ERR022075 # Illumina short reads
fi


