import matplotlib.pyplot as plt
from collections import Counter
import sys
import statistics  


if len(sys.argv) < 2:
    print("usage python3 utils/analyze_alignment <path to sl_raw_align.txt>")
    exit()




# # for i in [100000,200000]:  # Add more iteration values if needed
file_name = sys.argv[1]

# Read the file and extract the second column values
with open(file_name, "r") as file:
    F = [line.split(" ") for line in file]

C = Counter([f[1] for f in F])
C = [c for c in C.values() ]

mu = statistics.mean(C)

plt.hist(C,bins = 100)
plt.xlim(left = 0, right = mu*3)
plt.title(file_name)

plt.savefig("imgs/alignment_distrubution")



