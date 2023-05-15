from libgraph.basic_graph import *
from math import ceil

class CompactGraph:

    def __init__(self, old_graph = None):

        self.old_graph = old_graph
        self.new_graph = Graph()
        self.first_to_node = {}
        self.last_to_node = {}
        self.old_to_new = defaultdict(list) # maps old node id to new node ids
        self.new_to_old = defaultdict(set) # matches new node id to old node ids 
        self.node_size = -1

    def get_compact_graph(self, **kwargs):
        """
        Wrapper for different compacting functions
        
        Input:
            **argument --- a keyword dictionary. The first keyword must be "mode = <str>" specifying which mode of compaction is done.
            Other keywords are specified in different modes below.

        Current supporting modes:
        1. mode = branching
            - merges all 1-in-1-out nodes. No other keywords or arguments need to be provided
        """
        if "mode" not in kwargs:
            print("get_compact_graph: this function is not properly called. Must specify mode!")
            print("\t - branching")
            print("\t - empty_nodes, argument = aln_dict")
            print("\t - small_variations, argument = threshold")
            exit()

        mode = kwargs["mode"]

        if mode == "branching":
#            print("Merging non-branching nodes", flush=True)
            self.get_compact_graph_branching()
            self.old_to_new_path()
            self.validate_new_seq()


# ============================================================================================
#
#                            COMPACT MODE = REMOVE NON-BRANCHING NODES
#
# ============================================================================================

# ------
# Some design choices:
# Using the original logic on merging nodes to merge edges may induce parallel edges that have different labels.
# so instead of merging the entire branching path that starts and ends with the node with > 1 in/out degree, 
# we restrict the path to only start and end with 1-in-1-out nodes.

    def concat_edge_label(self, node_list) -> str:
        """
            Merge a list of nodes into new node The only function that updates old_to_new and new_to_old mapping.
            Note that the input is node but updates the edge.
            Basically contracts all the edges between the first and the last node.
            Args:
                node_list --- input list of nodes to be merged, node id in old_graph
            Return:
                concatenated edge label
        """
        s = ""
        weight = 0

        node1 = -1
        node2 = -1
        # print(node_list)
        if node_list[0] in self.old_to_new:
            node1 = self.old_to_new[node_list[0]][0]
        else:
            node1 = self.new_graph.add_node()
        
        if node_list[-1] in self.old_to_new:
            node2 = self.old_to_new[node_list[-1]][0]
        else:
            node2 = self.new_graph.add_node()
        
        for idx in range(len(node_list)-1):
            n1 = node_list[idx]
            n2 = node_list[idx + 1]
            s += self.old_graph.edge_labels[(n1, n2)]
            weight = self.old_graph.adj_list[n1][n2]
            self.old_to_new[n1] = [node1]
            self.new_to_old[node1].add(n1)

        self.old_to_new[node_list[-1]] = [node2]
        self.new_to_old[node2].add(node_list[-1])

        # print(node_list, s)

        self.new_graph.add_edge(node1, node2, s, count=weight,forbid_dup=True)
        
        return s

    def add_compact_node(self, non_branch_path):
        """
            Merge a non-branching path as a new compact node in the new graph.
        """
        self.first_to_node[non_branch_path[0]] = self.num_nodes
        self.last_to_node[non_branch_path[-1]] = self.num_nodes
        self.new_graph.add_node(self.concat_node_label(non_branch_path,self.num_nodes), self.num_nodes)
        self.num_nodes += 1
    
    def get_compact_graph_branching(self):

        source, sink = self.old_graph.source_in, self.old_graph.sink_out

        for v in self.old_graph.nodes:
            # skip one-in-one-out nodes. They will be taken care of by creating non-branching paths
            if self.old_graph.indeg[v] != 1 or self.old_graph.outdeg[v] != 1 or v==source or v==sink:
                for neighbor in self.old_graph.adj_list[v]:
                    self.concat_edge_label([v, neighbor])
                for neighbor in self.old_graph.adj_list_reverse[v]:
                    self.concat_edge_label([neighbor, v])

                non_branch_path = []                
                for u in self.old_graph.adj_list[v]:

                    while self.old_graph.indeg[u] == 1 and self.old_graph.outdeg[u] == 1 and u!=source and u!=sink:
                        non_branch_path.append(u)
                        u = list(self.old_graph.adj_list[u].keys())[-1]                    

                    if self.old_graph.indeg[u] == 1 and u!=sink:
                        non_branch_path.append(u)
                    if len(non_branch_path) > 1:
                        self.concat_edge_label(non_branch_path)

                    non_branch_path = []

    def old_to_new_path(self):
        """Convert paths in old graph to paths in new graph
        """
        try:
            new_paths = defaultdict(list)
            new_weights = defaultdict()
            for pathid in self.old_graph.paths:
                path = self.old_graph.paths[pathid]
                curr_new = self.old_to_new[path[0]][0]
                new_paths[pathid].append(curr_new)
                new_weights[pathid] = self.old_graph.weights[pathid]
                for i in range(len(self.old_graph.paths[pathid])):
                    if path[i] in self.new_to_old[curr_new]:
                        continue
                    curr_new = self.old_to_new[path[i]][0]
                    new_paths[pathid].append(curr_new)
            self.new_graph.paths = new_paths
            self.new_graph.weights = new_weights
        except Exception as e:
            print("old_to_new_path",e)

    def validate_new_seq(self):
        """validate the new graph contains the same sequences as the old graph
        """
        try:
            
            for seqid, seq in zip(self.old_graph.seqids, self.old_graph.seqs):
                path = self.new_graph.paths[seqid[0]]
                s = ""
                for i in range(len(path) - 2):
                    s += self.new_graph.edge_labels[(path[i], path[i+1])]
                assert(s==seq)
        except AssertionError:
            print("sequences not the same")
        except AttributeError as e:
            print("validate_new_seq", e)

        # for u in self.old_graph.adj_list:
        #     for v in self.old_graph.adj_list[u]:
        #         if u not in self.last_to_node or v not in self.first_to_node:
        #             continue
        #         n1 = self.last_to_node[u]
        #         n2 = self.first_to_node[v]
        #         self.new_graph.adj_list[n1][n2] += self.old_graph.adj_list[u][v]
        #         self.new_graph.adj_list_reverse[n2][n1] += self.old_graph.adj_list_reverse[v][u]
        #         self.new_graph.update_degree(n1, n2)
