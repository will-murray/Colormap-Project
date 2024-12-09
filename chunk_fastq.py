import sys
import os

folder = sys.argv[1]
fastq_file_1 = sys.argv[2]
fastq_file_2 = sys.argv[3]
chunk_size = int(sys.argv[4])  
MAX_CHUNKS = int(sys.argv[5])

def write_fastq_to_chunk_folder(fastq_filename):
    lines_in_chunk = 0
    chunk_num = 0
    chunk_lines = [] 


    with open(f"{folder}/{fastq_filename}", "r") as file:
        for line in file:
            if chunk_num > MAX_CHUNKS:
                return
            
            chunk_lines.append(line)  
            lines_in_chunk += 1

            if lines_in_chunk == chunk_size: 
                chunk_filename = f"{folder}/chunks/{fastq_filename}_{chunk_num}.fastq"
                with open(chunk_filename, "w") as chunk_file:
                    chunk_file.writelines(chunk_lines)  # Write the lines to the chunk file
                # print(f"Written chunk {chunk_num} to {chunk_filename}")

                chunk_num += 1
                lines_in_chunk = 0
                chunk_lines = []  

        if chunk_lines:
            chunk_filename = f"{folder}/chunks/{fastq_filename}_{chunk_num}.fastq"
            with open(chunk_filename, "w") as chunk_file:
                chunk_file.writelines(chunk_lines)
            print(f"Written final chunk {chunk_num} to {chunk_filename}")

write_fastq_to_chunk_folder(fastq_file_1)
write_fastq_to_chunk_folder(fastq_file_2)
