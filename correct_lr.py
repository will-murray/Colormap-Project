import sys
import re
import numpy as np
import time
import networkx as nx


if len(sys.argv) == 1: #test mode
    freq_file = "test_data/freq.txt"
    data_file = "test_data/sl_raw_align.txt"
    lr_file = "test_data/pac.fasta"
else:
    freq_file = sys.argv[1] #file containing the number of short reads mapped to each long read
    data_file = sys.argv[2] #file containing |short read name| long read name | alignment left end| alignment right end|
    lr_file = sys.argv[3] #long read file


def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def get_lr_dict(fname):
    D = {}
    with open(fname,"r") as file:
        contents = file.read()
        R = contents.split(">")
        for r in R[1:]:
           k,v = r.split('\n')[:2]
           D[k]=v
    
    return D


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


    



def build_long_read_graph(D, verbose = False):
    """
    builds the graph as an adjanceny matrix for a particular long read
    Input D: A numpy 2d array whose rows are the short reads mapped to the long read with name lr_name
    """

    lr_name = D[0,1]
    if lr_name == "*":
        return None
    
    G = nx.Graph()
    G.add_nodes_from([i[0] for i in D])

    print("-------------\n",lr_name) if verbose else 0

    candidates = 0
    found = 0
    # brute force (for now)
    for i in D:
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
                    w = levenshteinDistance(j[4][B[1]:], LR[lr_name][int(j[2]) + B[1] : int(j[2]) + B[1] + len(j[4][B[1]:])])
                    G.add_edge(i[0], j[0],weight= w)
                               
                    if verbose:
                        print(f"***********\n{i[0]} <-> {j[0]}")
                        print(f"{i[4][:A[0]].lower()}{i[4][A[0]:A[1]]} <-> {j[4][B[0]:B[1]]}{j[4][B[1]:].lower()}\n")
                        print(f"overlapping interval :{I}, interval length = {I[1]- I[0]}")
                        print(f"edge weight (edit distance): {w}")
                        print(f"read <-> ref | {j[4][B[1]:]} <-> { LR[lr_name][int(j[2]) + B[1] : int(j[2]) + B[1] + len(j[4][B[1]:])] }\n")

                    
                    found += 1

    
    return G

    

def build_graphs(fname):
    """
    returns a list of graphs, one for each long read
    """

    #init the data sets
    data = np.loadtxt(fname, dtype=str)
    
    
    #group the short read data based on the which long read its mapped to
    partitions = [data[start:end] for start,end in B]
    i = 0
    for D in partitions:
        print(f"{(i / len(partitions)) * 100}%")
        G = build_long_read_graph(D, verbose=False)
        print(G)
        

        i+=1

        if i == 4:
            exit()


    
            
B = compute_break_points(freq_file)
LR = get_lr_dict(lr_file)
MIN_OVERLAP = 10

t0 = time.time()
build_graphs(data_file)
print(f"runtime = {time.time() - t0}")