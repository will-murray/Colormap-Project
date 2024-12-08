# Snakefile to process input files and generate results.txt

# Output file

#           parameters          #

folder = "ecoli_new"
n_short_reads = 2000
n_long_reads = 100

#################################

n_short_reads *=4
n_long_reads *=2


output_file = f"{folder}/results.txt"


rule all:
    input:
        output_file

# Rule to process the FASTQ and FASTA files
rule align_short_to_long:
    output:
        alignment = f"{folder}/sl_align.sam",
        long_reads = f"{folder}/lr.fasta"
    
    shell:
        """
        awk 'NR <= {n_short_reads} {{print}}' ecoli_new/SRR13921543_1_sub.fastq > {folder}/sr1.fastq
        awk 'NR <= {n_short_reads} {{print}}' ecoli_new/SRR13921543_2_sub.fastq > {folder}/sr2.fastq
        awk 'NR <= {n_long_reads} {{print}}' ecoli_new/SRR10971019_sub.fasta    >  {folder}/lr.fasta 

        bwa index {folder}/lr.fasta
        bwa mem -t 4 -aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 10 -w 40 -r1 -D 0 -y 20 -L 30,30 -T 2.5 {folder}/lr.fasta {folder}/sr1.fastq {folder}/sr2.fastq > {output.alignment}



        """


rule format_short_long_alignment:
    input:
        alignment = f"{folder}/sl_align.sam"
    output:
        raw_alignment = f"{folder}/sl_raw_align.txt"
    shell:
        """

        grep -v "^@" {input.alignment} | #remove filter for non header lines from the same file
        awk '{{print $1 "." NR%2, $3, $4, $4 + length($10), $10}}' |
        sort -k2,2 -k3,3 -k4,4 | awk '$5 != "*" {{print}}' > {output.raw_alignment}
        
        """

rule correct_long_reads:
    input:
        raw_alignment = f"{folder}/sl_raw_align.txt",
        long_reads = f"{folder}/lr.fasta"
        
    output:
        output_file
    shell:
        """
    
        # g++ colormap.cpp -o colormap
        ./colormap {input.long_reads} {input.raw_alignment}
        echo "placeholder" > {output_file}
        """
    


