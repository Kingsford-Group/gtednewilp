'''
This grpah library is different in that the labels are on the edges instead of nodes. 
Additionally, the edges are stored in a multiset (the same edge could be stored multiple times)
'''

import sys
from collections import defaultdict, Counter
from Bio import SeqIO
import networkx as nx
# import vg_pb2
import stream
import pdb

class Graph:
    """
        A general genome graph class.
        self.nodes --- set of node indices
        self.edge_labels --- dictionary that maps (node1,node2) to edge label
        self.idx_to_edge --- dictionary that maps edge idx to pair of nodes
        self.edge_to_idx --- dictionary that maps an edge to a list of idxs
        self.adj_list --- dictionary of dictionary that maps each node to its successors and the count of the edge
        self.adj_list_reverse --- dictionary of dictionary that maps each node to its predecessors
        self.in/outdeg --- dict that maps node to its in/out degrees
    """

    def __init__(self, edges=None, fname=None, get_degrees=True):
        """Constructor that takes different graph informations

        Args:
            edges (dict, optional): Mappiong from (node1, node2) to a tuple (string, count). Defaults to None.
            fname (str, optional): gfa or vg file that describes a graph object. Defaults to None.
            get_degrees (bool, optional): Populates in/out deg. Defaults to True.
        """

        self.nodes = set()
        self.edge_labels = defaultdict()
        self.idx_to_edge = defaultdict()
        self.edge_to_idx = defaultdict(set)
        self.adj_list = defaultdict(Counter)
        self.adj_list_reverse = defaultdict(Counter)
        self.indeg = Counter()
        self.outdeg = Counter()
        self.indeg_weighted = Counter()
        self.outdeg_weighted = Counter()
        self.total_edges = 0

        self.source_in = -1
        self.sink_out = -1

        if edges is not None:
            for edge in edges:
                self.edge_labels[edge] = edges[edge][0]
                edge_count = edges[edge][1]
                self.adj_list[edge[0]][edge[1]] += edge_count
                self.adj_list[edge[1]][edge[0]] += edge_count

                self.indeg[edge[1]] += 1
                self.outdeg[edge[0]] += 1
                self.indeg_weighted[edge[1]] += edge_count
                self.outdeg_weighted[edge[0]] += edge_count

                for _ in range(edge_count):
                    self.idx_to_edge[self.total_edges] = edge
                    self.edge_to_idx[edge].add(self.total_edges)
                    self.total_edges += 1
                self.nodes.add(edge[0])
                self.nodes.add(edge[1])

        # needs to shift labels from nodes to edges
        if fname != None:
            suffix = fname.split(".")[-1]
            if suffix == "gfa":
                self.from_gfa_multi(fname) 
            elif suffix == "vg":
                self.from_vg(fname)

        # get degrees
        if get_degrees:
            self.indeg = defaultdict(int)
            self.outdeg = defaultdict(int)

            self.indeg_weighted = defaultdict(int)
            self.outdeg_weighted = defaultdict(int)

            self.get_degrees()

    def clear(self):
        self.nodes = set()
        self.edge_labels = defaultdict()
        self.idx_to_edge = defaultdict()
        self.edge_to_idx = defaultdict(set)
        self.adj_list = defaultdict(Counter)
        self.adj_list_reverse = defaultdict(Counter)
        self.indeg = Counter()
        self.outdeg = Counter()
        self.indeg_weighted = Counter()
        self.outdeg = Counter()
        self.total_edges = 0
        
        self.source_in = -1
        self.sink_out = -1


    def update_degree(self, node1, node2, count=1, new_edge=True):

        if new_edge:
            self.indeg[node2] += 1
            self.outdeg[node1] += 1

        self.indeg_weighted[node2] += count
        self.outdeg_weighted[node1] += count
        

    def get_reverse_adj(self):
        self.adj_list_reverse = defaultdict(Counter)
        for u in self.adj_list:
            for v in self.adj_list[u]:
                self.adj_list_reverse[v][u] = self.adj_list[u][v]


    def get_degrees(self):

        self.indeg = defaultdict(int)
        self.outdeg = defaultdict(int)

        self.indeg_weighted = defaultdict(int)
        self.outdeg_weighted = defaultdict(int)

        for node1 in self.adj_list:
            for node2 in self.adj_list[node1]:
                self.update_degree(node1, node2, count=self.adj_list[node1][node2])
        return self.indeg, self.outdeg


    def get_edge_seq(self, edge_id=-1, edge=None):
        """Get sequence of an edge by either edge idx or by edge in the form of (node1, node2)

        Args:
            edge_id (int, optional): edge idx. Defaults to 0.
            edge (tuple, optional): (node1, node2). Defaults to None.

        Returns:
            sequence: str
        """
        assert edge_id != -1 and edge != None, "get_edge_seq: needs one parameter: idx or edge"
        assert edge_id == -1 and edge == None, "get_edge_seq: needs one parameter: idx or edge"

        if edge != None:
            return self.edge_labels[edge]
        else:
            return self.edge_labels[self.idx_to_edge[edge_id]]


    def add_node(self):
        """
        Add node to graph

        Return: node id
        """
        node_id = len(self.nodes)
        self.nodes.add(node_id)
        return node_id


    def add_edge(self, node1, node2, seq="", count=1, forbid_dup = False):
        """Add edge to graph

        Args:
            node1 (int): _description_
            node2 (int): _description_
            seq (str, optional): edge sequence. In the case of duplicate edges, dont set this seq. Defaults to "".
            count (int, optional): number of copies of the added edge. Defaults to 1.
            forbid_dup (bool): whether duplicated edges are allowed
        """
        if node1 not in self.nodes or node2 not in self.nodes:
            print("Fail to add edge (%i,%i)--- nodes do not exist" % (node1, node2))
            return

        if forbid_dup and (node1, node2) in self.edge_to_idx:
            return
        
        if (node1, node2) in self.edge_to_idx:
            assert len(seq) == 0 or seq == self.edge_labels[(node1, node2)], "sequence mismatch!"
            self.adj_list[node1][node2] += count
            self.adj_list_reverse[node2][node1] += count

            for _ in range(count):
                self.edge_to_idx[(node1, node2)].add(self.total_edges)
                self.idx_to_edge[self.total_edges] = (node1, node2)
                self.total_edges += 1

            self.update_degree(node1, node2, count, new_edge=False)
        else:
            assert len(seq) > 0, "must specify a sequence"
            self.adj_list[node1][node2] += count
            self.adj_list_reverse[node2][node1] += count

            self.edge_labels[(node1, node2)] = seq

            for _ in range(count):
                self.edge_to_idx[(node1, node2)].add(self.total_edges)
                self.idx_to_edge[self.total_edges] = (node1, node2)
                self.total_edges += 1

            self.update_degree(node1, node2, count)


    def from_gfa(self, fname, skip=False):
        """_summary_

        Args:
            fname (_type_): _description_
            skip (bool, optional): Whether skip non ATCGst nodes. Defaults to False.
        """
        self.clear()
        labeled_nodes = defaultdict()
        for line in open(fname):
            l = line.rstrip("\n").split("\t")
            if "S" in l[0]:
                # skip nodes with only N or other unknown characters in it
                if ("A" in l[2] or "C" in l[2] or "T" in l[2] or "G" in l[2] or "s" in l[2] or "t" in l[2]) or (not skip):
                    labeled_nodes[int(l[1])] = l[2]
                    self.nodes.add(int(l[1]))
            if "L" in l[0]:
                n1 = int(l[1])
                n2 = int(l[3])
                if n1 in self.nodes and n2 in self.nodes:
                    self.add_edge(n1, n2, labeled_nodes[n2], 1)

    def from_gfa_multi(self, fname):
        """Read from gfa where the 6th field of edge line is the multiplicity

        Args:
            fname (str): path to gfa file
        """
        self.clear()
        self.paths = defaultdict()
        for line in open(fname):
            l = line.rstrip("\n").split("\t")
            if "S" in l[0]:
                self.nodes.add(int(l[1]))
            if "L" in l[0]:
                n1 = int(l[1])
                n2 = int(l[3])
                multiplicity = int(l[5])
                label = l[6]
                self.add_edge(n1, n2, label, multiplicity)
            if "P" in l[0]:
                pathid = l[1]
                path = [int(i) for i in l[2].split(",")]
                self.paths[pathid] = path


    def to_gfa_multi(self, fname):
        """
        Write to gfa where the 6th field of edge line is the multiplicity

        Args:
            fname (str): path to gfa file
        """

        with open(fname, "w") as ofile:
            ofile.write("\t".join(["H", "VN:Z:1.0"]))

            for n in self.nodes:
                node_seq = ""
                ofile.write("\n")
                ofile.write("\t".join(["S", str(n), node_seq]))

            for u in self.adj_list:
                for v in self.adj_list[u]:
                    ofile.write("\n")
                    ofile.write("\t".join(["L", str(u), "+", str(v), "+", str(self.adj_list[u][v]), self.edge_labels[(u,v)]]))


    def verify_eulerian(self):
        self.get_degrees()
        for v in self.indeg_weighted:
            if not self.indeg_weighted[v] == self.outdeg_weighted[v]:
                print(v)

    # def from_vg(self, fname):
    #     self.clear()

    #     # read protobuf format of vg
    #     vg_list = []
    #     with stream.open(fname, "rb") as istream:
    #         for data in istream:
    #             try:
    #                 vg = vg_pb2.Graph()
    #                 vg.ParseFromString(data)
    #                 vg_list.append(vg)
    #             except Exception:
    #                 continue

    #     for vg in vg_list:
    #         for node in vg.node:
    #             if "A" in node.sequence or "C" in node.sequence or "T" in node.sequence or "G" in node.sequence or "s" in node.sequence or "t" in node.sequence:
    #                 self.nodes[node.id] = node.sequence

    #     for vg in vg_list:
    #         for edge in vg.edge:
    #             n1 = getattr(edge, "from")
    #             n2 = int(edge.to)
    #             if n1 in self.nodes and n2 in self.nodes:
    #                 self.adj_list[n1][n2] += 1

    # def to_gfa(self, fname):
    #     with open(fname, "w") as ofile:
    #         ofile.write("\t".join(["H", "VN:Z:1.0"]) + "\n")

    #         for n in self.nodes:
    #             strings = self.nodes[n]
    #             if strings == "$":
    #                 strings = "s"
    #             if strings == "#":
    #                 strings = "t"
    #             ofile.write("\t".join(["S", str(n), strings]) + "\n")

    #         for u in self.adj_list:
    #             for v in self.adj_list[u]:
    #                 for _ in range(self.adj_list[u][v]):
    #                     ofile.write("\t".join(["L", str(u), "+", str(v), "+", "0M"]) + "\n")

    def get_nx_graph(self, flow=False):
        graph = nx.DiGraph()

        for u in self.adj_list:
            for v in self.adj_list[u]:

                if flow:
                    if self.nodes[u] == "t" and self.nodes[v] == "s":
                        continue
                    graph.add_edge(u, v, capacity=self.adj_list[u][v])
                else:
                    graph.add_edge(u, v, capacity=self.adj_list[u][v])
        return graph


    def get_nx_multigraph(self, flow=False):
        graph = nx.MultiDiGraph()

        for edge in self.edge_to_idx:
            for _ in range(len(self.edge_to_idx[edge])):
                graph.add_edge(edge[0], edge[1])

        return graph

    def vg_graph_add_source_sink(self, print_msg = False):
        # only use for graphs that don't have source or sink nodes but are DAGs

        sources = [v for v in self.nodes if self.indeg[v] == 0]
        sinks = [v for v in self.nodes if self.outdeg[v] == 0]

        if print_msg:
            print("Old sources", sources)
            print("Old sinks", sinks)

        # add source
        source = self.add_node("$")
        for s in sources:
            self.add_edge(source, s)

        # add sink
        sink = self.add_node("#")
        for s in sinks:
            self.add_edge(s, sink)

        return self
    
    def is_source(self, node):
        return self.nodes[node] == "$" or self.nodes[node] == "s"

    def is_sink(self, node):
        return self.nodes[node] == "#" or self.nodes[node] == "t"

    def get_source_sink(self):
        """Get source and sink node of the graph
            - If the degrees are not all even, there must exist a source and a sink with odd degrees
            - If the degrees are all even, look for the edge with "#" label that connects from sink to source
            - Otherwise, cannot determine if the graph has a source and a sink node annotated
        Returns:
            tuple(source, sink): a tuple of source and sink nodes. Returns (-1, -1) if cannot determine if the grpah has a source and sink node
        """
        if self.source_in != -1 and self.sink_out != -1:
            return self.source_in, self.sink_out
        
        # check if degrees are all even
        for n in self.nodes:
            if n not in self.indeg_weighted:
                self.source_in = n
            if n not in self.outdeg_weighted:
                self.sink_out = n
        if self.source_in != -1 and self.sink_out != -1:
            return self.source_in, self.sink_out
        
        sink_edge = [edge for edge in self.edge_labels if self.edge_labels[edge] == "#"]
        if len(sink_edge) == 1:
            return sink_edge[0][1], sink_edge[0][0]
        
        return -1, -1

def copy_adj_list(adj_list):
    new_adj_list = defaultdict(Counter)
    for u in adj_list:
        for v in adj_list[u]:
            new_adj_list[u][v] = adj_list[u][v]
    return new_adj_list
