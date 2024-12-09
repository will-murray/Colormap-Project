# This script downloads short and long read data for Escherichia_coli_K12_MG1655
#   
#   We are using data from the bio sample SAMN02604091, see here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN02604091&o=acc_s%3Aa



SRA_ID_SHORT="SRR13921546"
SRA_ID_LONG="SRR10971019"
OUTPUT_DIR="./"
SPLIT_FILES=true



####
echo "extracting reference genome"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/425/475/GCF_904425475.1_MG1655/GCF_904425475.1_MG1655_genomic.fna.gz
gunzip GCF_904425475.1_MG1655_genomic.fna.gz
cp GCF_904425475.1_MG1655_genomic.fna ref.fasta
rm GCF_904425475.1_MG1655_genomic.fna


# ###
# echo "extracting ${SRA_ID_LONG} (long) with fastq-dump..."
# fastq-dump ${SRA_ID_LONG} 
# echo "converting result to .fasta format..."
# awk '{if(NR%4==1) {print ">" substr($0, 2)} else if(NR%4==2) {print $0}}' ${OUTPUT_DIR}${SRA_ID_LONG}.fastq > ${OUTPUT_DIR}${SRA_ID_LONG}.fasta
# rm "${OUTPUT_DIR}${SRA_ID_LONG}.fastq"

# ###
# echo "Extracting ${SRA_ID_SHORT} (short reads) with fasterq-dump..."
# fastq-dump --split-files "${SRA_ID_SHORT}"



