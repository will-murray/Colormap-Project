import sys
import re
import numpy as np
import time
import networkx as nx
import graphviz as gv


if len(sys.argv) == 1: #test mode
    freq_file = "test_data/freq.txt"
    data_file = "test_data/sl_raw_align.txt"
    lr_file = "test_data/pac.fasta"
else:
    freq_file = sys.argv[1] #file containing the number of short reads mapped to each long read
    data_file = sys.argv[2] #file containing |short read name| long read name | alignment left end| alignment right end|
    lr_file = sys.argv[3] #long read file

def visualize_nx_graph(G : nx.Graph, name = "default"):
    """
    converts a nx.Graph object into a graphviz dot object, and produces a jpg titled with name
    """
    A = gv.Digraph(name)

    for n,data in G.nodes(data = True):
        if G.degree(n) > 0:
            lab = str(data.get('label',n))+ " | " + str(data.get('value', n)) # Use 'label' attribute if available, otherwise use node ID
            A.node(str(n), label=lab)

    for e in G.edges(data=True):
        A.edge(e[0], e[1],label=str(e[2]['weight']))
    

        

    A.render(filename=name, format="png", cleanup=True)
    


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

def build_long_read_graph(D, verbose = False) -> nx.Graph:
    """
    builds the graph as an adjanceny matrix for a particular long read
    Input D: A numpy 2d array whose rows are the short reads mapped to the long read with name lr_name
    """

    lr_name = D[0,1]
    if lr_name == "*":
        return nx.Graph()

    G = nx.Graph()
    for i in D:
        G.add_node(i[0], value = i[2])

    print("-------------\n",lr_name) if verbose else 0

    candidates = 0
    found = 0
    # brute force (for now)
    for idx,i in enumerate(D):
        for j in D[idx:]:
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
                    G.add_edge(i[0], j[0],weight= w, I = (int(i[2]),I[0]) )
                    # edge weight is the defined in the paper
                    # I stores:
                    #   1. the left endpoint of read i,
                    #   2. the start of the overlap between i and j
                    #   3. the right endpoint of read j
                    # I is used to construct a sequence of short reads which replaces a section of the long read
                    if verbose:
                        print(f"***********\n{i[0]} <-> {j[0]}")
                        print(f"{i[4][:A[0]].lower()}{i[4][A[0]:A[1]]} <-> {j[4][B[0]:B[1]]}{j[4][B[1]:].lower()}\n")
                        print(f"overlapping interval :{I}, interval length = {I[1]- I[0]}")
                        print(f"edge weight (edit distance): {w}")
                        print(f"read <-> ref | {j[4][B[1]:]} <-> { LR[lr_name][int(j[2]) + B[1] : int(j[2]) + B[1] + len(j[4][B[1]:])] }\n")

                    
                    found += 1

    
    
    return G

    
def SP_correction(G : nx.Graph,D, lr_name):
    """
    Given a networkx graph and the name of the long read:
    1. compute the shortest s-d path of each connected component
        s = the leftmost mapped node
        d = the rightmost mapped node
    2. replace base pairs in the long read with the reads that form shortest paths
    """
    if len(G.edges) > 20:
            visualize_nx_graph(G)
            for comp in list(nx.connected_components(G)):
                if len(comp) > 1:
                    C = [(node, int(G.nodes[node]['value'])) for node in comp]
                    source = min(C, key=lambda x: x[1])[0]
                    dest = max(C, key=lambda x: x[1])[0]
                    s_val = G.nodes[source].get('value')
                    d_val = G.nodes[dest].get('value')

                    if source != dest:
                        SP = nx.dijkstra_path(G, source, dest)

                        s = ""  # String to hold the corrected long read formed by SP
                        O = []  # Overlaps between sequences

                        # Compute overlaps based on edges in SP
                        for idx in range(len(SP) - 1):
                            I = G.get_edge_data(SP[idx], SP[idx + 1])['I']
                            O.append(I[1] - I[0])


                        

                        
                        
                        
                        print(f"{source} : {dest}")
                        print(f"{len(s)} | [{s_val}, {int(d_val)}]")                                    
                                            
            
                                
                                
                        print("-------")

                        
                        
            exit()



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
        SP_correction(G,D, lr_name=D[0,1])

        i+=1


    
            
B = compute_break_points(freq_file)
LR = get_lr_dict(lr_file)
MIN_OVERLAP = 5

t0 = time.time()
build_graphs(data_file)
print(f"runtime = {time.time() - t0}")