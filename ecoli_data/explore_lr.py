import re
import matplotlib.pyplot as plt

D = []
with open("lr_len_data.txt", "r") as file:
    for line in file:
        match = int( re.search("(?<=\=)\d+", line).group(0) )
        D.append(match)

print(f"min: {min(D)}, max: {max(D)}, average:{sum(D)/ len(D)}")