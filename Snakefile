
folder = "ecoli"
short_1 = "SRR13921546_1.fastq"
short_2 = "SRR13921546_2.fastq"
long = "SRR10971019.fasta"
ref = "ref.fasta"

# folder = "test_data"
# short_1 = "ill_1.fastq"
# short_2 = "ill_2.fastq"
# long = "pac.fasta"


test_name = "1" #no singleton corrections
correct_singletons = "no"
deg = 1


short_reads_per_chunk = 10
n_long_reads = 10
max_chunks = 2



#################################
if test_name != "":
    test_name = f"_{test_name}"
else:
    test_name = "_default"
    
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

        # if [ {deg} -eq 1 ]; then
        #     echo "simulating errors in the long reads"
        #     python3 degrade.py {folder}/lr.fasta
        #     cp {folder}/lr_deg.fasta {folder}/lr.fasta
        #     rm {folder}/lr_deg.fasta
        # fi
        
        # Step 2: Index the long reads FASTA file
        bwa index {folder}/lr.fasta

        # Step 3: Create the chunks for the short reads
        mkdir -p {folder}/chunks
        python3 utils/chunk_fastq.py {folder} {short_1} {short_2} {short_reads_per_chunk} {max_chunks}

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
            bwa mem -t 8 -aY -A 5 -B 11 -O 2,1 -E 4,3 -k 8 -W 16 -w 40 -r1 -D 0 -y 20 -L 30,30 -T 2.5 \
            {folder}/lr.fasta "$s1" "$s2" |
            grep -v "^@" | 
            awk '{{print $1 "." NR%2, $3, $4 - 1, $4 + length($10) -1, $10}}' |
            # awk '{{print $1 "." NR%2, $3, $4 - 1, $4 + length($10) -1, $10}}' |
            awk '$5 != "*" {{print}}' >> {output.raw_alignment} 

        done

        
        sort -k2,2 -k3,3 -k4,4 {output.raw_alignment} -o {output.raw_alignment}
        """


rule correct_long_reads:
    input:
        raw_alignment = f"{folder}/sl_raw_align.txt",
        long_reads = f"{folder}/lr.fasta"
        

    output:
        long_reads_corr = f"{folder}/lr_corr{test_name}.fasta"
    shell:
        """

        g++ colormap.cpp -o colormap
        ./colormap {input.long_reads} {input.raw_alignment} {correct_singletons}

        #add the length of each long read to the corrected long reads file
        awk '/^>/ {{if (seq) print length(seq); seq=""; header=$0}} !/^>/ {{seq=seq$0}} END {{if (seq) print length(seq)}}' {folder}/lr_corr.fasta |
        awk 'NR==FNR {{len[NR]=$1; next}} /^>/ {{print $0 " length=" len[++count]}} !/^>/ {{print}}' - {folder}/lr_corr.fasta > {folder}/tmp.fasta && mv {folder}/tmp.fasta {folder}/lr_corr.fasta

        cp {folder}/lr_corr.fasta {output.long_reads_corr} 
        rm {folder}/lr_corr.fasta

        echo -e "finished {test_name}" > {output_file}



        """


rule align_to_reference:
    input: 
        long_reads_corr = f"{folder}/lr_corr{test_name}.fasta",
        long_reads = f"{folder}/lr.fasta"

    output:
        output_file

    shell:
        """
            python3 utils/filt.py {input.long_reads_corr} {input.long_reads}
            blasr --header --bestn 1 {folder}/lr_filt.fasta {folder}/{ref} > {folder}/og.bam
            blasr --header --bestn 1 {input.long_reads_corr} {folder}/{ref} > {folder}/corr.bam

        

            python3 utils/analyze_out.py {folder}/og.bam    {input.long_reads} > {output_file}
            python3 utils/analyze_out.py {folder}/corr.bam  {input.long_reads_corr} >> {output_file}

            echo -e "finished {test_name}" > {output_file}


            
        """
