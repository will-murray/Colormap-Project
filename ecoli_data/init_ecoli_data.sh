# This script downloads the required files pertaining to Escherichia_coli_K12_MG1655, this includes:
# - reference genome : https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3/





SRA_ID_SHORT="SRR13921543"
SRA_ID_LONG="SRR10971019"
OUTPUT_DIR="./"
SPLIT_FILES=true


# echo "extracting ${SRA_ID_LONG} (long) with fastq-dump..."
# fastq-dump ${SRA_ID_LONG} 
# echo "converting result to .fasta format..."
# awk '{if(NR%4==1) {print ">" substr($0, 2)} else if(NR%4==2) {print $0}}' ${OUTPUT_DIR}${SRA_ID_LONG}.fastq > ${OUTPUT_DIR}${SRA_ID_LONG}.fasta
# rm "${OUTPUT_DIR}${SRA_ID_LONG}.fastq"


# echo "Extracting ${SRA_ID_SHORT} (short reads) with fasterq-dump..."
#     fastq-dump --split-files "${SRA_ID_SHORT}"



echo "Extracting the reference genome..."
    wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_904425475.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" > ref.fasta