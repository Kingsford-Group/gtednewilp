from basic_graph import *
import pickle, os
from copy import deepcopy
from compact_graph import *
import pdb

class DeBruijnGraph(Graph):
    '''
    De Bruijn graph construction methods. 
    DBG defined here is the edge-centric graph where labels are on the graphs.
    Each kmer is stored as a sequence of k edges where each edge has one charater in its label.
    
    Construction algorithm:
        - Create a node for each k-1 mer in the graph
        - create an edge between each adjacent pair of k-1 mers
        - it is possible that there is only one node that has imbalanced weights (|in - out| = 1)
        - Do DFS on the created graph:
            - for each not visited node:
                - if node is a source: create an edge for all letters in label(node)[:-1]
                - for all outgoing edges of node: assign label(node)[-1] to edge
                    - push neighbor in visit stack
    '''

    def __init__(self, seqids=None, seqs=None, k=None):

        # self.nodes = dict({0:"$",1:"#"})      # key: node_id, value: label
        # self.adj_list = defaultdict(dict)     # key: node_id, value: dict(n1_id: weight1, n2_id, weight2))
        # self.edge_to_idx = dict()             # key: (u,v), value: idx
        # self.idx_to_edge = dict()             # key: idx, value: (u, v, weight)
        # self.paths = [[0] for i in range(len(seqs))]
        # self.seqs = seqs
        # self.seqid = seqid

        super().__init__(get_degrees=False)
        
        self.source_in = 0
        self.sink_out = 1
        self.nodes = set([self.source_in, self.sink_out])

        self.dbg_adj_list = defaultdict(Counter)

        self.paths = defaultdict(list)
        self.seqs  = seqs
        self.seqids = seqids
        self.weights = dict()
        self.k = k
        self.kmers = dict()

        # self.init_source_sink()
        self.construct_kmers(seqids, seqs,k)
        # self.add_edge(self.sink_out, self.source_in, seq="#", count=self.outdeg_weighted[self.source_in])

    # def init_source_sink(self):
    #     self.nodes = set([0,1,2,3])
    #     self.add_edge(0,1,seq="s")
    #     self.add_edge(2,3,seq="t")


    def validate_paths(self):
        
        for seqid, seq in zip(self.seqids, self.seqs):
            path = self.paths[seqid[0]]
            s = ""
            for i in range(len(path) - 2):
                s += self.edge_labels[(path[i],path[i+1])]
            assert(s == seq)

    def add_edges(self, node_list, weight):
        for i in range(len(node_list)-1):
            edge_idx = len(self.idx_to_edge)
            u = node_list[i]
            v = node_list[i+1]
            edge = (u,v,weight)

            if v not in self.adj_list[u]:
                self.adj_list[u][v] = weight
                self.edge_to_idx[(u,v)] = edge_idx
                self.idx_to_edge[edge_idx] = edge
            else:
                self.adj_list[u][v] += weight
                edge_idx = self.edge_to_idx[(u,v)]
                self.idx_to_edge[edge_idx] = (u,v,self.adj_list[u][v])


    def construct_kmers(self, seqids, seqs, k):

        self.kmers = set()

        self.kmer_indeg = Counter()
        self.kmer_outdeg = Counter()
        self.kmer_indeg_weighted = Counter()
        self.kmer_outdeg_weighted = Counter()


        for idx, seq in zip(seqids, seqs):
            weight = int(idx[1])
            next_kmer = ""
            for i in range(len(seq) - k + 1):
                curr_kmer = seq[i : i+k-1]
                next_kmer = seq[i+1 : i+k]

                self.kmers.add(curr_kmer)
                self.kmers.add(next_kmer)
                
                self.dbg_adj_list[curr_kmer][next_kmer] += weight
                self.kmer_indeg[next_kmer] += 1
                self.kmer_outdeg[curr_kmer] += 1
                self.kmer_indeg_weighted[next_kmer] += weight
                self.kmer_outdeg_weighted[curr_kmer] += weight

        i = 0
        with open("dbg.test.gfa", "w") as ofile:
            node_dict = dict()
            for kmer in self.kmers:
                ofile.write("S" +"\t" + str(i) + "\t" + kmer +"\n")
                node_dict[kmer] = i
                i += 1
            for kmer in self.dbg_adj_list:
                for kmer2 in self.dbg_adj_list[kmer]:
                    ofile.write("L" + "\t" + str(node_dict[kmer]) + "\t+" + "\t" + str(node_dict[kmer2]) + "\t+\t" +str(self.dbg_adj_list[kmer][kmer2])+"M" + "\n")
            for s in self.seqs:
                ofile.write(s +"\n")
        #for kmer in self.kmer_indeg_weighted:
            #if kmer != seqs[0][-k+1:]:
                #if self.kmer_indeg_weighted[kmer] != self.kmer_outdeg_weighted[kmer]:
                    #print(kmer, self.kmer_indeg_weighted[kmer],self.kmer_outdeg_weighted[kmer], self.dbg_adj_list[kmer])
        

        # sources = []
        # for kmer in self.kmers:
        #     if self.kmer_outdeg[kmer] - self.kmer_indeg[kmer] == 1:
        #         sources.append((kmer, self.kmer_outdeg_weighted[kmer] - self.kmer_indeg_weighted[kmer]))
        # print(self.kmer_outdeg)
        # print(self.kmer_indeg)
        # print(self.kmer_outdeg_weighted)
        # print(self.kmer_indeg_weighted)
        self.convert_graph()

        # for idx, seq in zip(seqid, seq):

        #     weight = int(idx[1])
        #     for i in range(len(seq) - k + 1):
        #         curr_kmer = seq[i : i+k-1]
        #         next_kmer = seq[i+1 : i+k]
         
                    
    def convert_graph(self):
        """
        Convert node-based k-1-mer graphs into an edge-labeled Euleiran graph
        """
        self.kmer_to_node = defaultdict()
        self.node_to_kmer = defaultdict()

        visited = set()
        for idx, seq in zip(self.seqids, self.seqs):
            source = seq[:self.k-1]
            path = self.split_node(source, int(idx[1]))
            # self.add_edge(self.source_in, first_node, seq="s", count=int(idx[1]))
            self.dfs_convert(source, visited)
            self.add_edge(self.kmer_to_node[seq[-(self.k-1):]], self.sink_out, seq="t", count=int(idx[1]))
            self.weights[idx[0]] = int(idx[1])
            self.paths[idx[0]] = path
            for i in range(len(seq)-self.k+2):
                curr_kmer = seq[i : i+self.k-1]
                self.paths[idx[0]].append(self.kmer_to_node[curr_kmer])
            self.paths[idx[0]].append(self.sink_out)


    def dfs_convert(self, kmer, visited):
        """
        Do DFS in the graph. push labels on the nodes to adjacent edges.

        kmer: current kmer
        node: last node on kmer's path
        """

        stack = [kmer]
        while len(stack) > 0:
            kmer = stack.pop()
            if kmer not in visited:
                visited.add(kmer)
                node = self.kmer_to_node[kmer]
                for next_kmer in self.dbg_adj_list[kmer]:
                    next_node = -1
                    if next_kmer not in self.kmer_to_node:
                        next_node = len(self.nodes)
                        self.nodes.add(next_node)
                        self.kmer_to_node[next_kmer] = next_node
                        self.node_to_kmer[next_node] = next_kmer
                    else:
                        next_node = self.kmer_to_node[next_kmer]
                    self.add_edge(node, next_node, seq=next_kmer[-1], count=self.dbg_adj_list[kmer][next_kmer]) 
                    # self.dfs_convert(next_kmer, visited)
                    stack.append(next_kmer)

                # if len(self.dbg_adj_list[kmer]) == 0:
                #     self.add_edge(node, self.sink_out, seq="t", count=self.kmer_indeg_weighted[kmer])


    def split_node(self, kmer, weight):
        """ kmer is split into a sequence of labeled edges
        Return a path of nodes
        """
        path = []
        n1 = self.source_in
        # self.nodes.add(n1)

        prev_node = n1
        for i in range(len(kmer)-1):
            path.append(prev_node)
            n2 = len(self.nodes)
            self.nodes.add(n2)
            self.add_edge(prev_node, n2, seq=kmer[i], count=weight)
            prev_node = n2
        path.append(prev_node)
        n2 = -1
        if kmer not in self.kmer_to_node:
            n2 = len(self.nodes)
            self.nodes.add(n2)
            self.kmer_to_node[kmer] = n2
            self.node_to_kmer[n2] = kmer
        else:
            n2 = self.kmer_to_node[kmer]
        
        self.add_edge(prev_node, n2, seq=kmer[-1], count=weight)
        return path

    def to_gfa_multi(self, fname):
        super().to_gfa_multi(fname)
        with open(fname, "a") as ofile:
            for path_idx in self.paths:
                ofile.write("\n")
                ofile.write("\t".join(["P"] + [str(self.weights[path_idx])] + [",".join([str(i) for i in self.paths[path_idx]])]))


def read_fasta(fname, get_seq=False):
    string_set = [("_".join(seq.id.split("_")[:-1]), seq.id.split("_")[-1]) for seq in SeqIO.parse(fname, "fasta")]
    seqs = []
    if get_seq:
        seqs = [str(seq.seq) for seq in SeqIO.parse(fname, "fasta")]
        return string_set, seqs
    return string_set

def test():
    seq = ["GACGG", "CGGCG"]
    seqid = [["1", "3"], ["2", "2"]]

    graph = DeBruijnGraph(seqid, seq, 3)
    print(graph.kmers)
    print(graph.nodes)
    print(graph.idx_to_edge)
    print(graph.edge_to_idx)
    print(graph.adj_list)
    print(graph.paths)
    print(graph.edge_labels)
    graph.to_gfa_multi("test.gfa")
    graph.verify_eulerian()



if __name__ == "__main__":

    # test()
    
    if len(sys.argv) < 4:
        print("USAGE: python3.6 %s <fasta_dir> <k-mer size> <graph output dir>" % sys.argv[0])
        exit()

    fasta_dir = sys.argv[1]
    k = int(sys.argv[2])
    target_dir = sys.argv[3]


    graph = Graph()

    for i in range(9):
        graph.add_node()
    graph.source_in = 0
    graph.source_out = 1

    graph.add_edge(0,2,"A",1)
    graph.add_edge(2,3,"B",1)
    graph.add_edge(3,4,"C",1)
    graph.add_edge(4,1,"D",1)
    graph.add_edge(0,5,"E",1)
    graph.add_edge(5,6,"F",1)
    graph.add_edge(6,7,"G",1)
    graph.add_edge(7,1,"H",1)
    graph.add_edge(0,8,"I",1)
    graph.add_edge(8,3,"J",1)

    prefix = "test"


    files = [f for f in os.listdir(fasta_dir) if ".fasta" in f]

    f = files[0]
    prefix=f.split(".")[0]

    seq_id, seqs = read_fasta(fasta_dir + "/" + f, True)
    seq_id = [("0","1")]
    graph = DeBruijnGraph(seq_id, seqs, k)
    graph.verify_eulerian()

    compact_graph = CompactGraph(old_graph=graph)
    compact_graph.get_compact_graph(mode="branching")

    # # print(compact_graph.new_to_old)
    # # print(compact_graph.old_to_new)


    graph_ofile = target_dir + "/" + prefix + "." + str(k) + ".dbg" 
    graph.to_gfa_multi(graph_ofile+".gfa")

    compact_graph.new_graph.to_gfa_multi(graph_ofile+".compact.gfa")

    # for f in os.listdir(fasta_dir):
    #     if ".fasta" in f:
    #         prefix=f.split(".")[0]
    #         print(f, prefix)
    #         seq_id, seqs = read_fasta(fasta_dir + "/" + f, True)
    #         seq_id = [("0","1")]
    #         graph = DeBruijnGraph(seq_id, seqs, k)
    #         # print(len(graph.nodes), len(graph.idx_to_edge))
    #         # with open(target_dir+"/"+prefix+"."+str(k)+".dbg.graph", "wb") as handle:
    #         #     pickle.dump(graph, handle)
    #         graph_ofile = target_dir + "/" + prefix + "." + str(k) + ".dbg.gfa"
    #         graph.to_gfa_multi(graph_ofile)
             
