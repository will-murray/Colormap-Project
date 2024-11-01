import sys
import re
import numpy as np
if len(sys.argv) == 1: #test mode
    freq_file = "test_data/freq.txt"
    data_file = "test_data/sl_raw_align.txt"
else:
    freq_file = sys.argv[1] #file containing the number of short reads mapped to each long read
    data_file = sys.argv[2] #file containing |short read name| long read name | alignment left end| alignment right end|



def build_freq_table(fname):
    T = {}
    with open(fname,"r") as file:
        for line in file:
            k,v = line.split(" ")
            v = re.sub("\n","",v)
            T.update({k:int(v)})

    return T

def build_long_read_graph(D, lr_name):
    if lr_name == "*":
        print("BRUH")

    

def build_graphs(fname):
    """
    returns a list of graphs, one for each long read
    """

    #init the data sets
    data = np.loadtxt(fname, dtype=str)

    current_lr_num = data[0,1]
    start_idx = 0
    end_idx = T[current_lr_num] #first row w/ data for the kth long read
    N = data.shape[0]

    while(end_idx  < N):
        # print(f"long read {current_lr_num} : {start_idx}->{end_idx}")
        build_long_read_graph(data[start_idx:end_idx])
        start_idx = end_idx
        current_lr_num = data[start_idx,1]
        end_idx += T[current_lr_num]

    # print(f"N = {N}, end_idx = {end_idx}")

    
            
T = build_freq_table(freq_file)
build_graphs(data_file)
