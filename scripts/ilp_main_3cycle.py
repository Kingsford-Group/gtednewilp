'''
Simulate 3cycle graphs and run CCTED and GTED ILPs on them.

Select 100 3cycle graphs with GTED - CCTED ranging from 1 to 5.
'''

import sys, os
sys.path.append("/home/yutongq/gtednewilp/libgraph")

from scripts.ilp_writers import *
from scripts.alignment_graph import *
from libgraph.basic_graph import *
from libgraph.DeBruijnGraph import *

from tqdm import tqdm

import os
import pandas as pd
from collections import defaultdict
import gurobipy as gp

from numpy.random import poisson
from numpy.random import normal
import numpy.random as random 
import matplotlib.pyplot as plt

import time

gp.setParam("OutputFlag", 0)
gp.setParam("LogToConsole", 0)
out_dict = defaultdict(list)

num_cycles = 3
len_mean = 5
overlap_mean = 1

def simulate_seq(length):
    s = ""
    if length > 0:
        
        for i in range(length):
            s += np.random.choice(["A","C","G","T"])
    return s

def simulate_3cycle_graph(overlap_space_idx, overlaps, lengths, seqs, overlap_strings, middle_string):
    graph = Graph()

    # add edges for middle cycle
    node_idx = 0
    graph.add_node()
    for i, s in enumerate(middle_string):
        if i == len(middle_string) - 1:
            graph.add_edge(i, 0, seq=s)
        else:
            graph.add_node()
            graph.add_edge(i, i+1, seq=s)

    # add edges for other cycles
    node_idx = overlaps[0]
    for i, seq in enumerate(seqs):
        init_node = node_idx
        for j, s in enumerate(seq[:lengths[i] - overlaps[i]]):
            if j == lengths[i] - overlaps[i] - 1:
                graph.add_edge(node_idx, init_node - overlaps[i], s)
                node_idx = init_node - overlaps[i]
            else:
                graph.add_node()

                graph.add_edge(node_idx, len(graph.nodes)-1, seq=s)
                node_idx = len(graph.nodes)-1
        for ii in range(overlaps[i]):
            graph.add_edge(node_idx, node_idx + 1, graph.edge_labels[(node_idx, node_idx + 1)])
            node_idx += 1

        if i < len(overlap_space_idx1):
            init_node = node_idx + overlap_space_idx[i] + overlaps[i+1]
            node_idx = init_node
    return graph


for fii in tqdm(range(100)):
    
    while True:

        try:
            overlap_space_idx1,overlap_space_idx2, overlaps, lengths, seqs, overlap_strings,middle_string1=0,0,0,0,0,0,""

            while(middle_string1 == ""):
                try:
                    lengths = [int(normal(len_mean, 2)) for i in range(num_cycles)]
                    lengths[-1] = 15
                    overlaps = [abs(int(normal(overlap_mean, 1))) for i in range(num_cycles-1)]
                    print(lengths, overlaps)

                    seqs = [simulate_seq(lengths[i]) for i in range(num_cycles-1)]

                    overlap_space_idx1 = sorted([np.random.choice(lengths[i]-overlaps[i]) for i in range(num_cycles-2)])
                    overlap_space_idx2 = sorted([np.random.choice(lengths[i]-overlaps[i]) for i in range(num_cycles-2)])
                    overlap_strings = [seqs[i][len(seqs[i]) - overlaps[i]:] for i in range(num_cycles-1)]

                    middle_string1 = ""
                    for i in range(num_cycles-2):
                        middle_string1 += overlap_strings[i] + simulate_seq(overlap_space_idx1[i])
                    middle_string1 += overlap_strings[-1] + simulate_seq(max(lengths[-1]-sum(overlaps)-sum(overlap_space_idx1), 3))
                except Exception as e:
                    print(e)

            g1 = simulate_mickymouse_graph(overlap_space_idx1, overlaps, lengths, seqs, overlap_strings,middle_string1)
            g2 = simulate_mickymouse_graph(overlap_space_idx2, overlaps, lengths, seqs, overlap_strings,middle_string1)

            ag = Alignment_Graph(g1, g2)


            t0 = time.time()
            model1 = CCTED_ilp(ag, g1, g2, "test.lp", True)
            ccted = model1.getObjective().getValue()

            ccted_time = time.time() - t0
            
            t0 = time.time()
            model2 = GTED_ilp(ag, g1, g2, "test.lp", True)
            gted = model2.getObjective().getValue()
            
            if  gted - ccted == (fii// 20 + 1):
                print(model2.getObjective().getValue())
                out_dict["gted2"].append(gted)
                out_dict["time_gted_compact"].append(time.time() - t0)

                out_dict["ccted"].append(ccted)
                out_dict["time_ccted"].append(ccted_time)
                break
                
        except Exception as e:
            print(e)

    g1.to_gfa_multi("data/3cycle_1_%i.gfa" % fii)
    g2.to_gfa_multi("data/3cycle_2_%i.gfa" % fii)
    ag = Alignment_Graph(g1, g2)

    t0 = time.time()
    num_it = 1
    model1 = CCTED_ilp(ag, g1, g2, "test.lp", True)
    
    while True:
        model1,subgraph1 = update_scc_constr(ag,model1)
        if subgraph1 == -1:
            out_dict["iteration"].append(num_it)
            out_dict["gted"].append(model1.getObjective().getValue())
            break
        model1.reset()
        model1.optimize()
        if model1.status != GRB.OPTIMAL:
            print(model1.status)
            out_dict["iteration"].append(num_it)
            out_dict["gted"].append(-1)
            break
        if num_it % 50 == 0:
            print(num_it, model1.getObjective().getValue(), flush=True)
        if time.time() - t0 > 1200:
            out_dict["iteration"].append(-num_it)
            out_dict["gted"].append(model1.getObjective().getValue())
            break
        num_it += 1
    out_dict["time_gted_it"].append(time.time() - t0)
    print( model1.getObjective().getValue())
    
    if fii % 10 == 9:
        df = pd.DataFrame(out_dict)
        df.to_csv("3cycle_output.csv", index=False)