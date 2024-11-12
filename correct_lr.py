import sys
import re
import numpy as np
import time
import statistics as stat

if len(sys.argv) == 1: #test mode
    freq_file = "test_data/freq.txt"
    data_file = "test_data/sl_raw_align.txt"
else:
    freq_file = sys.argv[1] #file containing the number of short reads mapped to each long read
    data_file = sys.argv[2] #file containing |short read name| long read name | alignment left end| alignment right end|



def compute_break_points(fname):
    T = []
    with open(fname,"r") as file:
        for line in file:
            k,v = line.split(" ")
            v = int(re.sub("\n","",v)) #cast to int, remove newline
            T.append([k,v])
    
    T.sort(key= lambda x: x[0])
    B = []
    start_idx = 0
    for t in T:
        B.append([start_idx, start_idx + t[1]])
        start_idx += t[1]


    return B


    



def build_long_read_graph(D):
    """
    builds the graph as an adjanceny matrix for a particular long read
    Input D: A numpy 2d array whose rows are the short reads mapped to the long read with name lr_name
    """

    lr_name = D[0,1]
    if lr_name == "*":
        return 0
    G = []
    candidates = 0
    found = 0
    # brute force (for now)
    for i in D:
        edges = []
        for j in D:
            self_edge =   i[0] == j[0]
            i_before_j =  int(i[2]) <= int(j[2]) and int(i[3]) < int(j[3])
            overlapping = int(j[2]) <= int(i[3]) - MIN_OVERLAP + 1

            if not self_edge and i_before_j and overlapping:
                candidates += 1
                I = [int(j[2]), int(i[3])] #overlapping interval
                A = [I[0]-int(i[2]) , I[1] -int(i[2])]
                B = [I[0]-int(j[2]) , I[1] -int(j[2])]

                if( i[4][A[0]:A[1]] == j[4][B[0]:B[1]]):
                    # # print(f"overlapping interval :{I}, interval length = {I[1]- I[0]}")
                    # print(i[4][A[0]:A[1]] == j[4][B[0]:B[1]])
                    found += 1
    if candidates == 0:
        return 0
    return found/candidates

    

def build_graphs(fname):
    """
    returns a list of graphs, one for each long read
    """

    #init the data sets
    data = np.loadtxt(fname, dtype=str)
    
    
    #group the short read data based on the which long read its mapped to
    partitions = [data[start:end] for start,end in B]
    S = []
    i = 0
    for p in partitions:
        print(f"{(i / len(partitions)) * 100}%")
        S.append(build_long_read_graph(p))
        i+=1

    print(f"mean : {stat.mean(S)}")
    print(f"stddev : {stat.stdev(S)}")
    print(f"min : {min(S)}")
    print(f"max : {max(S)}")
    return

    
            
B = compute_break_points(freq_file)
MIN_OVERLAP = 10

t0 = time.time()
build_graphs(data_file)
print(f"runtime = {time.time() - t0}")