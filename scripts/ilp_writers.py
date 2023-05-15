import gurobipy as gp
from gurobipy import GRB
from scripts.alignment_graph import *
from libgraph.basic_graph import *
from libgraph.compact_graph import *
from libgraph.DeBruijnGraph import *
import networkx as nx
import Levenshtein

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


def CCTED_ilp(aln_graph, graph_1, graph_2, ofilename, ilp=True):
    '''
    Objective: min cost(e) f(e)
    s.t. flow symmetry and flow coverage constraints
    ilp = True if computes ILP, else, all variables are continuous
    '''
    try:
        # create variables set
        E = [idx for idx in aln_graph.idx_to_edge]  # each variable corresponds to an edge in alignment graph

        # create a new model
        model = gp.Model("GTED")
        
        # creat variables
        x = 0
        if ilp:
            x = model.addVars(E, vtype=GRB.INTEGER, name="x")
        else:
            x = model.addVars(E, vtype=GRB.CONTINUOUS, name="x")
            
        # set objective
        model.setObjective(sum([aln_graph.idx_to_edge[idx][2] * x[idx] for idx in E]), GRB.MINIMIZE)

        # source and sink
        source = aln_graph.nodes[0]
        sink = aln_graph.nodes[1]
        model.addConstr(sum([x[i] for i in outgoing_edges(aln_graph, source)]) == 1)
        model.addConstr(sum([x[i] for i in incoming_edges(aln_graph, sink)]) == 1)

        # flow conservation/symmetry
        for node in aln_graph.nodes.values():
            if node != source and node != sink:
                model.addConstr(sum([x[i] for i in incoming_edges(aln_graph, node)]) - sum([x[i] for i in outgoing_edges(aln_graph, node)]) == 0)

        # flow coverage
        for e1 in aln_graph.project_1:
            edge1 = graph_1.idx_to_edge[e1]
            total_weight = graph_1.adj_list[edge1[0]][edge1[1]]
            model.addConstr(sum([x[i] for i in aln_graph.project_1[e1]]) == total_weight)
        for e2 in aln_graph.project_2:
            edge2 = graph_2.idx_to_edge[e2]
            total_weight = graph_2.adj_list[edge2[0]][edge2[1]]
            model.addConstr(sum([x[i] for i in aln_graph.project_2[e2]]) == total_weight) 

        model.optimize()

        if model.status == GRB.INTERRUPTED:
            print("!!! INTERRUPTED")
            return model

        if model.status == GRB.OPTIMAL:
            return model

        if model.status == GRB.INFEASIBLE:
            print("!!! INFEASIBLE MODEL?")
            model.write(ofilename)
            return model

        if model.status == GRB.TIME_LIMIT:
            print("!!! TIME LIMIT EXCEEDS")
            model.write(ofilename)
            return model
    
    except gp.GurobiError as e:
        sys.exit('!!! GurobiError: Error code ' + str(e.errno) + ': ' + str(e))    


def GTED_ilp(aln_graph, graph_1, graph_2, ofilename, ilp=True):
    '''
    Objective: min cost(e) f(e)
    s.t. flow symmetry and flow coverage constraints
    '''
    try:
        # create variables set
        E = [idx for idx in aln_graph.idx_to_edge]  # each variable corresponds to an edge in alignment graph

        # create a new model
        model = gp.Model("GTED")
        
        # creat variables
        x = 0
        d = 0
        y = 0
        if ilp:    
            x = model.addVars(E, vtype=GRB.INTEGER, name="x")
            d = model.addVars(E, vtype=GRB.INTEGER, name="d")
            y = model.addVars(E, vtype=GRB.BINARY, name="y")
        else:
            x = model.addVars(E, vtype=GRB.CONTINUOUS, name="x")
            d = model.addVars(E, vtype=GRB.CONTINUOUS, name="d")
            y = model.addVars(E, vtype=GRB.CONTINUOUS, name="y")

        # set objective
        model.setObjective(sum([aln_graph.idx_to_edge[idx][2] * x[idx] for idx in E]), GRB.MINIMIZE)

        # source and sink
        source = aln_graph.nodes[0]
        sink = aln_graph.nodes[1]
        model.addConstr(sum([x[i] for i in outgoing_edges(aln_graph, source)]) == 1)
        model.addConstr(sum([x[i] for i in incoming_edges(aln_graph, sink)]) == 1)

        # flow conservation/symmetry
        for node in aln_graph.nodes.values():
            if node != source and node != sink:
                model.addConstr(sum([x[i] for i in incoming_edges(aln_graph, node)]) - sum([x[i] for i in outgoing_edges(aln_graph, node)]) == 0)

        # flow coverage
        for e1 in aln_graph.project_1:
            edge1 = graph_1.idx_to_edge[e1]
            total_weight = graph_1.adj_list[edge1[0]][edge1[1]]
            model.addConstr(sum([x[i] for i in aln_graph.project_1[e1]]) == total_weight)
        for e2 in aln_graph.project_2:
            edge2 = graph_2.idx_to_edge[e2]
            total_weight = graph_2.adj_list[edge2[0]][edge2[1]]
            model.addConstr(sum([x[i] for i in aln_graph.project_2[e2]]) == total_weight) 

        # ordering constraint
        model.addConstr(sum([d[i] for i in outgoing_edges(aln_graph, source)]) == 1)
        model.addConstr(sum([d[i] for i in incoming_edges(aln_graph, sink)]) == sum([x[i] for i in range(len(x))]))
        for node in aln_graph.nodes.values():
            if node != source and node != sink:
                model.addConstr(sum([d[i] for i in outgoing_edges(aln_graph, node)]) - sum([d[i] for i in incoming_edges(aln_graph, node)]) == sum([x[i] for i in outgoing_edges(aln_graph, node)]))

        # if x[i] = 0, then d[i] = 0
        for i in aln_graph.idx_to_edge:
            model.addConstr(-x[i] - len(aln_graph.idx_to_edge) * y[i] <= -1)
            model.addConstr(d[i] - len(aln_graph.idx_to_edge) * (1-y[i]) <= 0)
            model.addConstr(d[i] >= 0)
            model.addConstr(x[i]>=0)
            model.addConstr(y[i]>=0)
            model.addConstr(y[i]<=1)

        model.optimize()

        if model.status == GRB.INTERRUPTED:
            print("!!! INTERRUPTED")
            return model

        if model.status == GRB.OPTIMAL:
            return model

        if model.status == GRB.INFEASIBLE:
            print("!!! INFEASIBLE MODEL?")
            model.write(ofilename)
            return model

        if model.status == GRB.TIME_LIMIT:
            print("!!! TIME LIMIT EXCEEDS")
            model.write(ofilename)
            return model
    
    except gp.GurobiError as e:
        sys.exit('!!! GurobiError: Error code ' + str(e.errno) + ': ' + str(e))    


def construct_sol_graph(ag,model):
    '''
    Construct subgraph induced by solutions in a ILP model
    '''
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



