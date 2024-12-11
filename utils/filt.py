import sys

if len(sys.argv) != 3:
    print("usage:python3 filt.py corr.fasta og.fasta")
    exit()

corr_fname,og_fname = sys.argv[1], sys.argv[2]


with open(corr_fname, "r") as file:
    targets = [read.split(" ")[0] for read in file if read.startswith(">")]    


newfile = og_fname[:len(og_fname) - 6] + "_filt.fasta"
write_seq = False
with open(og_fname, "r") as file, open(newfile,"w") as nf:
    for line in file:
        if(line.startswith(">")):
            if(line.split(" ")[0] in targets):
                nf.write(line)
                write_seq = True
        elif write_seq == True:
            nf.write(line)
            write_seq = False
                
    