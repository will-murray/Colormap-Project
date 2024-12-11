import sys

if len(sys.argv) != 3:
    print("usage: python3 utils/add_missing_lr <origonal long reads>.fasta <corrected long reads>.fasta")

og_fname = sys.argv[1]
corr_fname = sys.argv[2]


OG = {}
full_keys = {}
with open(og_fname, "r") as file:
    key = ""
    full_key = ""
    val = ""
    for line in file:
        if line.startswith(">"):
            key = line.split(" ")[0]
            full_key = line
        else:
            val = line
            OG[key] = val
            full_keys[key] = full_key

with open(corr_fname, "r") as file:
    corrected_keys = [line.split(" ")[0] for line in file if line.startswith(">")]

missing_keys = [key for key in OG.keys() if key not in corrected_keys]

with open(corr_fname, "a") as file:
    for key in missing_keys:
        # print(full_keys[key])
        # print(OG[key])
        file.write(full_keys[key])
        file.write(OG[key])
