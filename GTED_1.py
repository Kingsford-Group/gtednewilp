import sys, os
#sys.path.append("/home/yutongq/gtednewilp/libgraph")

from basic_graph import *
from DeBruijnGraph import *
from FGTED_utils_updated import *
from GTED_main import *
import gurobipy as gp
from tqdm import tqdm
import Levenshtein
import time
import pandas as pd
import pdb
import gurobipy as gp
from tqdm import tqdm
import Levenshtein
import time
from Bio import Align
from numpy.random import poisson
from numpy.random import normal
import numpy.random as random 
import matplotlib.pyplot as plt

import networkx as nx


def select_seed(s,length,mut_rate,not_loc = -1):
    loc = int(not_loc)
    while loc in range(not_loc, not_loc + length):
        loc = random.choice(len(s)-length - 1)
#     print(loc)
    s = s[loc: loc+length]
    return s, int(loc)

def insert_mutation(s, mut_rate):
    new_s = ""
    for i in range(len(s)):
        if random.choice(10) < mut_rate:
            mut = random.choice(["A", "T", "C", "G", "-", "-"])
            if mut != "-":
                new_s += mut
        else:
            new_s += s[i]
    return new_s


def get_str_dist(graph1, graph2):
    s1 = ""
    path1 = graph1.paths['1']
    path2 = graph2.paths['1']

    for i in range(len(graph1.paths['1'])-1):
        s1 += graph1.edge_labels[graph1.paths['1'][i], graph1.paths['1'][i+1]]

    s2 = ""
    for i in range(len(graph2.paths['1'])-1):
        s2 += graph2.edge_labels[graph2.paths['1'][i], graph2.paths['1'][i+1]]

    return Levenshtein.distance(s1,s2), s1, s2

def construct_sol_graph(ag,model):
    # get all x variables
    xs = dict()
    for v in model.getVars():
        if round(v.x) > 0:
            if "x" in v.VarName:
                idx = int(v.VarName.lstrip("x[").rstrip("]"))
                xs[idx] = round(v.x)
    
    subgraph = nx.MultiDiGraph()
    for idx in xs:
        n1,n2,c,w,e_type = ag.idx_to_edge[idx]
        for _ in range(xs[idx]):
#             if n1 == n2:
#                 print(n1, n2)
            subgraph.add_edge(n1,n2)
            
    subgraph.add_edge(ag.nodes[1], ag.nodes[0])
    
    return subgraph, xs     


def update_scc_constr(ag,model):
    subgraph,xs = construct_sol_graph(ag,model)
    sccs = []
    if nx.number_strongly_connected_components(subgraph)==1:
        return model,-1
    for s in nx.strongly_connected_components(subgraph):
        if ag.nodes[1] not in s and ag.nodes[0] not in s:
            sccs.append(s)
    '''
    if sccs == []:
        return model, -1
    '''
    for scc in sccs:
        scc_edges = set()
        scc_out = []
        #scc_in = []
        total_weight = 0
        scc_size = 0
        for e in xs:
            n1, n2, c, w, e_type = ag.idx_to_edge[e]
            if n1 in scc and n2 in scc:
                scc_edges.add(e)
                scc_size +=1
        scc_size = int(scc_size)
        for e in ag.idx_to_edge:
            n1, n2, c, w, e_type = ag.idx_to_edge[e]
            if n1 in scc and e not in scc_edges:
                scc_out.append(e)
            if e in scc_edges:
                total_weight+=w
        total_weight = int(total_weight)
        x = [v for v in model.getVars() if "x" in v.VarName]
        y = model.addVar(vtype=GRB.BINARY)
        #         total_weight = len(scc)
        model.addConstr(sum(x[e] for e in scc_edges) >= scc_size *y)
        model.addConstr(sum(x[e] for e in scc_edges) - scc_size + 1 - total_weight * y <= 0)
        model.addConstr(sum([x[e] for e in scc_out]) >= y)
        #model.addConstr(sum([x[e] for e in scc_in]) >= y)
        model.update()
    return model, subgraph


'''
data_dir = "data/"
graph_files = [data_dir + f for f in os.listdir(data_dir) if "4.gfa" in f]

seg_len = 15
mut_rate = 0
num_copies = 3
nucs = []

length1 = int(normal(seg_len, 2))
length2 = int(normal(seg_len, 2))
# mut_rate = poisson(mut_rate)
# num_copies1 = poisson(num_copies)
num_copies1 = 2

i = 0
out_dict = defaultdict(list)

# for file1 in graph_files[1:5]:
j = 0
#     for file2 in graph_files[1:5]: 
file1 = graph_files[0]
file2 = graph_files[0]
prefix1 = file1.split("/")[-1].split(".")[0] 
prefix2 = file2.split("/")[-1].split(".")[0]

graph1 = Graph(fname = file1)
graph2 = Graph(fname = file2)

ed, s1, s2 = get_str_dist(graph1, graph2)

s1_seed1, loc = select_seed(s1, length1,5)
s1_seed2, loc = select_seed(s1, length1,5, loc)
s2_seed1, loc = select_seed(s2, length2,5)
s2_seed2, loc = select_seed(s2, length2,5, loc)

s1_seed1 = "BBBDDDBBB"

locs = [[2, 15], [9, 30]]

str_list = [s1,s2]
for i,s in enumerate(str_list):
    inserted_locs =set([0])
    inserted_locs2 = []
    new_seed = []
    for seed in [s1_seed1]:
        #, s1_seed2, s2_seed1, s2_seed2]:
        for _ in range(num_copies1):
            new_s = insert_mutation(seed, mut_rate)
#             loc = random.choice(100)
# #                 print(loc, inserted_locs[i])
#             while loc in inserted_locs:
#                 loc = int(random.choice(100))
#             for l in range(loc, loc+len(new_s)):
#                 inserted_locs.add(l)
#             inserted_locs2.append(loc)
            new_seed.append(new_s)

    total_increased = 0
    seed_loc = sorted([(i,j) for i,j in zip(new_seed, locs[i])], key=lambda x: x[1])
    for seed, loc in seed_loc:
        print(seed, loc)
        new_loc = loc + total_increased
        str_list[i] = str_list[i][:new_loc] + seed + str_list[i][new_loc:]
        total_increased += len(seed)
    print(seed_loc)
s11 = str_list[0][:50]
s22 = str_list[1][:50]
print(s11)
print(s22)
print(ed, Levenshtein.distance(s11, s22), len(s11), len(s22))

g1 = DeBruijnGraph(seqids = [(1,1)], seqs=[s11], k=4)
g2 = DeBruijnGraph(seqids = [(1,1)], seqs=[s22], k=4)

out_dict["s1"].append(s11)
out_dict["s2"].append(s22)

g1.to_gfa_multi(data_dir+prefix1 + "_mut_%i_%i" % (i,j)+"_1" + ".gfa")
g2.to_gfa_multi(data_dir+prefix2 + "_mut_%i_%i" % (i,j) +"_2" ".gfa")

#start = time.time()

ag = Alignment_Graph(g1, g2)

t0 = time.time()
ilp=True
model1 = CCTED_ilp(ag, g1, g2, "test.lp", ilp)
print("wall clock time: "+str(time.time()-t0))


while True:
    model1,subgraph1 = update_scc_constr(ag,model1)
    if subgraph1 == -1:
        break
    #nx.draw(subgraph1)
    model1.reset()
    model1.optimize()
'''