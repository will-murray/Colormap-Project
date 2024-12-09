# Snakefile to process input files and generate results.txt

# Output file

#           parameters          #
# awk '/^>/ {if (seq) print length(seq); seq=""; header=$0} !/^>/ {seq=seq$0} END {if (seq) print length(seq)}' pac_corr.fasta | awk 'NR==FNR {len[NR]=$1; next} /^>/ {print $0 "_length=" len[++count]} !/^>/ {print}' - pac_corr.fasta > output.fasta
#add this to the pipeline


folder = "ecoli"
short_1 = "SRR13921543_1.fastq"
short_2 = "SRR13921543_2.fastq"
long =    "SRR10971019_sub.fasta"


test_name = "nsc" #no singleton corrections
correct_singletons = "no"
short_reads_per_chunk = 5000
n_long_reads = 50000
max_chunks = 1000

#################################

short_reads_per_chunk *=4
n_long_reads *=2


og_alignment = f"{folder}/og.sam"
corr_alignment = f"{folder}/corr.sam"

output_file = f"{folder}/results.txt"


rule all:
    input:
        output_file

# Rule to process the FASTQ and FASTA files
rule align_short_to_long:
    output:
        raw_alignment = f"{folder}/sl_raw_align.txt",
        long_reads = f"{folder}/lr.fasta"
    shell:
        """
        awk 'NR <= {n_long_reads} {{print}}' {folder}/{long} > {folder}/lr.fasta 
        
        # Step 2: Index the long reads FASTA file
        bwa index {folder}/lr.fasta

        # Step 3: Create the chunks for the short reads
        mkdir -p {folder}/chunks
        python3 chunk_fastq.py {folder} {short_1} {short_2} {short_reads_per_chunk} {max_chunks}

        # Step 4: Calculate number of chunk files
        num_files=$(ls -l {folder}/chunks/ | grep -v '^d' | wc -l)
        num_files=$(((num_files - 1) / 2))

        # Step 5: Loop through the chunk files and run bwa mem for each pair
        for ((i=0; i<num_files; i++)); do
            # Reference the chunk files using $i
            s1="{folder}/chunks/{short_1}_$i.fastq"
            s2="{folder}/chunks/{short_2}_$i.fastq"
            
            echo "Processing: $s1, $s2"
            
            # Run bwa mem for each chunk of short reads
            bwa mem -v 1 -t 8 -aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 10 -w 40 -r1 -D 0 -y 20 -L 30,30 -T 2.5 \
            {folder}/lr.fasta "$s1" "$s2" |
            grep -v "^@" | 
            awk '{{print $1 "." NR%2, $3, $4, $4 + length($10), $10}}' |
            awk '$5 != "*" {{print}}' >> {output.raw_alignment} 

        done

        
        sort -k2,2 -k3,3 -k4,4 {output.raw_alignment} -o {output.raw_alignment}
        """



# rule format_short_long_alignment:
#     input:
#         alignment = f"{folder}/sl_align.sam"
#     output:
#         raw_alignment = f"{folder}/sl_raw_align.txt"
#     shell:
#         """

#         grep -v "^@" {input.alignment} | #remove filter for non header lines from the same file
#         awk '{{print $1 "." NR%2, $3, $4, $4 + length($10), $10}}' |
#         sort -k2,2 -k3,3 -k4,4' > {output.raw_alignment}
        
#         """

rule correct_long_reads:
    input:
        raw_alignment = f"{folder}/sl_raw_align.txt",
        long_reads = f"{folder}/lr.fasta"
        

    output:
        output_file
    shell:
        """
    
        g++ colormap.cpp -o colormap
        ./colormap {input.long_reads} {input.raw_alignment} {correct_singletons}

        #add the length of each long read to the corrected long reads file
        awk '/^>/ {{if (seq) print length(seq); seq=""; header=$0}} !/^>/ {{seq=seq$0}} END {{if (seq) print length(seq)}}' {folder}/lr_corr.fasta |
        awk 'NR==FNR {{len[NR]=$1; next}} /^>/ {{print $0 " length=" len[++count]}} !/^>/ {{print}}' - {folder}/lr_corr.fasta > {folder}/tmp.fasta && mv {folder}/tmp.fasta {folder}/lr_corr.fasta

        cp {folder}/lr_corr.fasta {folder}/lr_corr_{test_name}.fasta 
        rm {folder}/lr_corr.fasta

        echo -e "finished {test_name}" > {output_file}



        """




        
        
        
        
        
        
        
        
        
        # echo -e "og" > {output.og_alignment}
        # echo -e "corr" > {output.corr_alignment}