from GTED_1 import *
import os
import argparse
import pandas as pd
from collections import defaultdict
import gurobipy as gp
import random
import pdb

gp.setParam("OutputFlag", 0)
gp.setParam("LogToConsole", 0)

parser = argparse.ArgumentParser()
parser.add_argument('--file1', type=str)
parser.add_argument('--file2', type=str)
args = parser.parse_args()

data_dir = "data/"

fnames = [f for f in os.listdir(data_dir) if ".4." in f ]

out_dict = defaultdict(list)

random.shuffle(fnames)

#fnames = ["simulated262.4.gfa","simulated280.4.gfa"]
for i,f1 in enumerate(fnames):
    for f2 in fnames[i:]:

        g1 = Graph(fname = data_dir+f1)
        g2 = Graph(fname = data_dir+f2)

        prefix1=f1.split(".")[0]
        prefix2=f2.split(".")[0]
        print(prefix1, prefix2)

        out_dict["g1"].append(prefix1)
        out_dict["g2"].append(prefix2)

        ag = Alignment_Graph(g1, g2)

        t0 = time.time()
        ilp=True
        model1 = CCTED_ilp(ag, g1, g2, "test.lp", ilp)
        out_dict["ccted"].append(model1.getObjective().getValue())

        # print("wall clock time: "+str(time.time()-t0))
        ccted_time = time.time()-t0

        t0 = time.time()
        num_it = 1
        while True:
            #pdb.set_trace()
            model1,subgraph1 = update_scc_constr(ag,model1)
            if subgraph1 == -1:
                out_dict["iteration"].append(num_it)
                out_dict["gted"].append(model1.getObjective().getValue())
                break
            #get_opt_aln(subgraph1,ag,g1,g2)
            #nx.draw(subgraph1)
            model1.reset()
            model1.optimize()
            if model1.status != GRB.OPTIMAL:
                print(model1.status)
                out_dict["iteration"].append(num_it)
                out_dict["gted"].append(-1)
                break
            print(num_it, flush=True)
            num_it += 1
        # print("wall clock time: "+str(time.time()-t0))
        gted_time = time.time() - t0
        ed, s1, s2 = get_str_dist(g1,g2)
        out_dict["edit_distance"].append(ed)
        out_dict["gted_time"].append(gted_time+ccted_time)
        out_dict["ccted_time"].append(ccted_time)
        print(prefix1,prefix2,ed,out_dict["ccted"][-1], out_dict["gted"][-1], ccted_time,gted_time+ccted_time, num_it)


df = pd.DataFrame(out_dict)
df.to_csv("ilp_output_more_0502.csv", index=False)
