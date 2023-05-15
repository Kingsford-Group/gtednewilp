import numpy as np
from collections import defaultdict
from Bio import SeqIO
from collections import Counter
import pandas as pd
import sys
import time

rev_dict = {"A":"T", "C":"G","G":"C","T":"A"}
gap_penalty = 1
mismatch_penalty = 1
match_reward = 0
large_penalty = 1
  
def matched(label_1, label_2):
    '''
    Determine match score
    1) if characters: return match or mismatch
    2) free edges from both sink to both source, otherwise forbidden
    '''
    if label_1 == label_2:
        return 0
    else:
        return 1
    
def find_gap(u1,u2,v1,v2):

    '''
    Impose large gap penalty if any node involves a sink
    In other words, only allow edges to sink if it is a match edge
    '''
    return gap_penalty

class Alignment_Graph(object):
    
    def __init__(self, graph_1, graph_2):
        
        self.nodes = dict()                  # key: (u1, u2) value: idx
        self.project_1 = defaultdict(set)    # maps edge from graph 1 to its projection on alignment graph
        self.project_2 = defaultdict(set)    # maps edge from graph 2 to its projection on alignment graph
        self.edge_to_idx = defaultdict(list)            # key: (n1,n2,cost,weight) value: idx
        self.idx_to_edge = dict()            # key: idx value:(n1,n2,cost, weight)
        self.idx_to_node = dict()            # key: idx value:(u1, u2)
        
        self.incoming = defaultdict(list)    # key: node_id1, value: [node_id2] for edge (node2, node1)
        self.outgoing = defaultdict(list)    # key: node_id1, value: [node_id2] for edge (node1, node2)

        self.graph_1 = graph_1
        self.graph_2 = graph_2
        
        node_counter = 0
        # create aligned nodes
        for u1 in graph_1.nodes:
            for u2 in graph_2.nodes:
                self.nodes[(u1,u2)] = node_counter
                self.idx_to_node[node_counter] = (u1, u2)
                node_counter += 1

        edge_counter = 0
        # create aligned edges
        for u1 in graph_1.nodes:
            for v1 in graph_1.adj_list[u1]:
                for u2 in graph_2.nodes:
                    for v2 in graph_2.adj_list[u2]:
                        
                        # obtain edge in original graphs
                        edges1 = list(graph_1.edge_to_idx[(u1,v1)])
                        edges2 = list(graph_2.edge_to_idx[(u2,v2)])
                        weight1 = graph_1.adj_list[u1][v1]
                        weight2 = graph_2.adj_list[u2][v2]
                        label1 = graph_1.edge_labels[(u1,v1)]
                        label2 = graph_2.edge_labels[(u2,v2)]
                        
                        max_weight = max(weight1, weight2)

                        # horizontal edge
                        h_edge = (self.nodes[(u1,u2)], self.nodes[(v1,u2)], gap_penalty, weight1, "d")

                        # vertical edge
                        v_edge = (self.nodes[(u1,u2)], self.nodes[(u1,v2)], gap_penalty, weight2, "d")

                        # diagonal edge
                        cost = matched(label1,label2)
                        d_edge = (self.nodes[(u1,u2)], self.nodes[(v1,v2)], cost, max_weight, "d")

                        edge1 = edges1[0]
                        edge2 = edges2[0]
                        edge_counter = self.create_edges([h_edge,v_edge,d_edge], edge_counter, edge1, edge2)
                
        # add source and sink edges
        source = 0
        sink = 1
        self.nodes[source] = node_counter
        self.idx_to_node[node_counter] = source
        node_counter += 1

        self.nodes[sink] = node_counter
        self.idx_to_node[node_counter] = sink
        node_counter += 1

        # euleri = 0 if closed_trail = 2 if open_trail, exit if other number
        euler1,s1,t1 = isEulerian(graph_1)
        euler2,s2,t2 = isEulerian(graph_2)
        assert euler1 == 0 or euler1 == 2
        assert euler2 == 0 or euler2 == 2

        # if both closed Eulerian trail
        if euler1 == 0 and euler2 == 0:
            s1 = list(graph_1.nodes)[0]
            s2 = list(graph_2.nodes)[0]
            for node1 in graph_1.nodes:
                edge_counter = self.add_edge(source, (node1, s2), 0, 1, edge_counter, "ss")
                edge_counter = self.add_edge((node1,s2), sink, 0, 1, edge_counter, "ss")
            for node2 in graph_2.nodes:
                if node2 == s2: continue
                edge_counter = self.add_edge(source, (s1, node2), 0, 1, edge_counter, "ss")
                edge_counter = self.add_edge((s1,node2), sink, 0, 1, edge_counter, "ss")

        elif euler1 ==2 and euler1 == 2:
            edge_counter = self.add_edge(source, (s1,s2), 0, 1, edge_counter, "ss")
            edge_counter = self.add_edge((t1,t2), sink, 0, 1, edge_counter, "ss")
        else:
            if euler1 == 2 and euler2 == 0:
                for node in euler2.nodes:
                    edge_counter = self.add_edge(source, (s1, node),  0, 1, edge_counter, "ss")
                    edge_counter = self.add_edge((t1, node), sink, 0, 1, edge_counter, "ss")
            else:
                for node in euler1.nodes:
                    edge_counter = self.add_edge(source, (node, s2), 0, 1, edge_counter, "ss")
                    edge_counter = self.add_edge((node, t2), sink, 0, 1, edge_counter, "ss")

    def add_edge(self, n1, n2, cost, weight, edge_counter, e_type):
        nn1 = self.nodes[n1]
        nn2 = self.nodes[n2]
        edge = (nn1, nn2, cost, weight,e_type)

        self.idx_to_edge[edge_counter] = edge
        self.edge_to_idx[edge].append(edge_counter)
        self.incoming[edge[1]].append(edge_counter)
        self.outgoing[edge[0]].append(edge_counter)

        edge_counter += 1
        return edge_counter

    def check_all_edges(self):
        for idx in self.idx_to_edge:
            edge1,edge2,cost,weight = self.idx_to_edge[idx]
            u1,u2 = self.idx_to_node[edge1]
            v1,v2 = self.idx_to_node[edge2]
            label_1 = self.graph_1.edge_labels[(u1,v1)]
            label_2 = self.graph_2.edge_labels[(u2,v2)]

            if cost == 0 and (label_1 != label_2):
                print("!!", u1,v1,u2,v2, label_1, label_2)
                        
    def edge_to_string(self, edge_idx):
        edge1,edge2,cost,weight = self.idx_to_edge[edge_idx]
        u1,u2 = self.idx_to_node[edge1]
        v1,v2= self.idx_to_node[edge2]
        label_1 = self.graph_1.nodes[v1]
        label_2 = self.graph_2.nodes[v2]
        return "[(%i,%i),(%i,%i)],%i,%i,%s,%s" % (u1,u2,v1,v2,cost,weight,label_1,label_2), "{%i,%i}"%(u1,u2),"{%i,%i}"%(v1,v2),"{{%i,%i}->{%i,%i}, \"%s,%s,%i " % (u1,u2,v1,v2,label_1,label_2,cost), 

    def create_edges(self, edges, edge_counter, e1, e2):
        '''
        Updates the data structures with newly created AG edges
        Parameters: edges -- [h_edge, v_edge, d_edge] ~ must follow the same order!
                    edge_counter -- keep tracks of the number of edges
                    e1, e2 -- edges in the original edges
        '''
        for i,edge in enumerate(edges):
            
            # don't add edge if edge has large penalty
            #if edge[2] == large_penalty:
            #    continue

            # check if edge already exists
            if edge not in self.edge_to_idx:
                self.idx_to_edge[edge_counter] = edge
                self.edge_to_idx[edge].append(edge_counter)
                self.incoming[edge[1]].append(edge_counter)
                self.outgoing[edge[0]].append(edge_counter)

                edge_idx = edge_counter
                edge_counter += 1

                # update projection
                if i == 0: # horizontal
                    self.project_1[e1].add(edge_idx)
                if i == 1: # vertical
                    self.project_2[e2].add(edge_idx)
                if i == 2: # diagonal
                    self.project_1[e1].add(edge_idx)
                    self.project_2[e2].add(edge_idx)

        return edge_counter

def isEulerian(graph):
    '''
    Find the number of nodes that are not degree-balanced in graph
    Returns:
    - euler: number of nodes that are not degree-balanced
    - s: node where in < out
    - t: node where in > out
    '''
    euler = 0
    s = -1
    t = -1
    for node in graph.nodes:
        if node in graph.indeg_weighted and node in graph.outdeg_weighted:
            if graph.indeg_weighted[node] > graph.outdeg_weighted[node]:
                t = node
                euler += 1
            if graph.indeg_weighted[node] < graph.outdeg_weighted[node]:
                s = node
                euler += 1
        else:
            euler += 1
            if node in graph.indeg_weighted:
                t = node
            if node in graph.outdeg_weighted:
                s = node
    return euler, s, t

def incoming_edges(aln_graph, node):
    '''
    helper function for FGTED_lp
    outputs a list of edge id for edges coming into node
    '''
    return [n for n in aln_graph.incoming[node]]

def outgoing_edges(aln_graph, node):
    '''
    helper function for FGTED_lp
    outputs a list of edge id for edges coming out of node
    '''
    return [n for n in aln_graph.outgoing[node]]