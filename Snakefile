
# folder = "ecoli"
# short_1 = "SRR13921546_1.fastq"
# short_2 = "SRR13921546_2.fastq"
# long = "SRR10971019.fasta"
# ref = "ref.fasta"

# folder = "test_data"
# short_1 = "ill_1.fastq"
# short_2 = "ill_2.fastq"
# long = "pac.fasta"
# ref = None

folder = "yeast"
short_1 = "ERR10616375_1.fastq"
short_2 = "ERR10616375_2.fastq"
long = "SRR18210286.fasta"
ref = "ref.fasta"

test_name = "1" #no singleton corrections
correct_singletons = "no"
deg = 1


short_reads_per_chunk = 10000
n_long_reads = 1000
max_chunks = 5



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


        
        #Index the long reads FASTA file
        bwa index {folder}/lr.fasta

        #Create the chunks for the short reads
        mkdir -p {folder}/chunks
        python3 utils/chunk_fastq.py {folder} {short_1} {short_2} {short_reads_per_chunk} {max_chunks}

        #Calculate number of chunk files
        num_files=$(ls -l {folder}/chunks/ | grep -v '^d' | wc -l)
        num_files=$(((num_files - 1) / 2))

        #Loop through the chunk files and run bwa mem for each pair
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

        
        sort -k2,2 -k3,3n -k4,4n {output.raw_alignment} -o {output.raw_alignment}
        """


rule correct_long_reads:
    input:
        raw_alignment = f"{folder}/sl_raw_align.txt",
        long_reads = f"{folder}/lr.fasta"
        

    output:
        long_reads_corr = f"{folder}/lr_corr{test_name}.fasta"
    shell:
        """

        g++ fast_colormap.cpp -o fast_colormap
        ./fast_colormap {input.long_reads} {input.raw_alignment} {correct_singletons}

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
            if [ "{folder}" -eq "test_data" ]; then
                echo "No reference genome for this data; no report produced" >> {output_file}
                exit 1
            fi
            # If colormap.cpp didnt correct any base pairs in a long read L then L wont be written to input.long_reads_corr,
            # this script add these missing (unchanged) reads
            # this is done to make the testing consistent

            python3 utils/add_missing_lr.py {input.long_reads} {input.long_reads_corr}

            blasr --header --bestn 1 {input.long_reads} {folder}/{ref} > {folder}/og.bam
            blasr --header --bestn 1 {input.long_reads_corr} {folder}/{ref} > {folder}/corr.bam

            
            echo -e "short reads per chunk: {short_reads_per_chunk}" >> {output_file}
            echo -e "num chunks: {max_chunks}">> {output_file}
            echo -e "number of long reads corrected: {n_long_reads}">> {output_file}
            echo -e "correct singletons?: {correct_singletons}">> {output_file}

            python3 utils/analyze_out.py {folder}/og.bam    {input.long_reads} >> {output_file}
            python3 utils/analyze_out.py {folder}/corr.bam  {input.long_reads_corr} >> {output_file}


            
        """
