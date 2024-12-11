import sys
import random as rd

if(len(sys.argv) != 2):
    print("Usage: python3 degrade.py <long_reads>.fasta")
fasta_file = sys.argv[1]
B = ["A","C","T","G"]

F = []
with open(fasta_file, "r") as file:
    for line in file:
        line_mut = ""
        if line.startswith(">"):
            F.append(line)
        else:
            for bp in line:
                if rd.random() < 0.15:
                    bp = B[rd.randint(0,3)]
                line_mut += bp
            F.append(line_mut)

deg_file = fasta_file[:len(fasta_file) - 6] + "_deg.fasta"
with open(deg_file ,"x") as file:
    for line in F:
        file.write(line)
