
from scripts.ilp_writers import *
from scripts.alignment_graph import *
from libgraph.basic_graph import *
from libgraph.DeBruijnGraph import *

import os
import pandas as pd
from collections import defaultdict
import gurobipy as gp
import random

gp.setParam("OutputFlag", 0)
gp.setParam("LogToConsole", 0)

data_dir = "data/"

fnames = [f for f in os.listdir(data_dir) if "3cycle_1_" in f ]
out_dict = defaultdict(list)
fnames = sorted(fnames)

ilp=True   # change to False to run LPs

for i, f1 in enumerate(fnames):
    idx = f1.split(".")[0].split("_")[-1]
    f2 = "3cycle_2_" + idx + ".gfa"

    g1 = Graph(fname = data_dir + f1)
    g2 = Graph(fname = data_dir + f2)

    ag = Alignment_Graph(g1, g2)

    # ccted ilp
    t0 = time.time()
    model1 = CCTED_ilp(ag, g1, g2, "test.lp", ilp)
    out_dict["ccted"].append(model1.getObjective().getValue())
    ccted_time = time.time()-t0
    out_dict["ccted_time"].append(ccted_time)

    # gted compact ilp
    t0 = time.time()
    model2 = GTED_ilp(ag, g1, g2, "test.lp", True)
    gted = model2.getObjective().getValue()
    out_dict["gted_compact"].append(gted)
    out_dict["time_gted_compact"].append(time.time() - t0)

    # gted iterative ilp
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
        if time.time() - t0 > 1200:
            # use negative iteration number to indicate early termination
            out_dict["iteration"].append(-num_it)
            out_dict["gted"].append(model1.getObjective().getValue())
            break
        num_it += 1
    out_dict["time_gted_it"].append(time.time() - t0)

df = pd.DataFrame(out_dict)
df.to_csv("3cycle_output.csv", index=False)