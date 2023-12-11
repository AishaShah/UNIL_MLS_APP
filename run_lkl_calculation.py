import argparse
import numpy as np
import scipy
from scipy.linalg import expm, sinm, cosm
import math

class TreeNode:
    def __init__(self, node_name):
        self.seq_leaf = None
        self.children = []
        self.parent = None
        self.branch_length = None
        self.lkl_list = None
        self.name = node_name

    def __len__(self):
        return len(self.seq_leaf)

    def getSequence(self):
        return self.seq_leaf

    def getBranchLength(self):
        return self.branch_length

    def getChildrenNodes(self):
        return self.children

    def get_lkl_list(self):
        return self.lkl_list

    def initialize_lkl_list(self):
        if self.seq_leaf:
            bases = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1]}
            self.lkl_list = [bases[base] for base in self.seq_leaf]
        else:
            self.lkl_list = None

def add_nodes_recursively(node, parent, node_dict, bl=None, 
                          sequence=None, lkl_list=None):
    if node not in node_dict:
        new_node = TreeNode(node)
        node_dict[node] = new_node

    if parent is not None:
        if parent not in node_dict:
            add_nodes_recursively(node=parent, parent=None, node_dict=node_dict)
        
        # Now, parent should be in node_dict
        node_dict[parent].children.append(node_dict[node])
        node_dict[node].parent = node_dict[parent]
        node_dict[node].branch_length = bl
        node_dict[node].seq_leaf = sequence
        node_dict[node].lkl_list = []
        node_dict[node].initialize_lkl_list()

# Reading tree topology, branch lengths, and sequences
def read_input_files(tree_path, branch_lengths_path, msa_path):
    with open(tree_path, "r") as tree_file, open(branch_lengths_path, "r") as BL_file, open(msa_path, "r") as msa_file:
        line_num = 0
        bl = BL_file.readline().split(",")
        sequence = {}
        for line in msa_file:
            key, val = line.split()
            sequence[key] = val
        
        # Dictionary to store references to nodes by their names
        node_dict = {}

        for parent_child in tree_file:
            tmp_list = parent_child.split(",")
            parent_node_name = tmp_list[0]
            node_name = tmp_list[1].rstrip()
            branch_len = bl[line_num]

            if node_name in sequence.keys():
                leaf_seq = sequence[node_name]
            else:
                leaf_seq = None

            add_nodes_recursively(node=node_name,
                                  parent=parent_node_name,
                                  node_dict=node_dict,
                                  bl=float(branch_len),
                                  sequence=leaf_seq)

            line_num += 1

    return node_dict

def calculate_lkl_anc(childA, childB, Q, lkl_dict):
    childA_lkl = lkl_dict[childA.name] # we will get a list of lists with liklihood for each nuclleotide at each position
    childB_lkl = lkl_dict[childB.name]
    anc_AB_lkl = []

    # calculating the liklihood of all 4 nucleotides at each position for the parent
    # using all positions in both children nodes, calculate ancetral nucleotide liklihoods
    for pos_list in range(len(childA_lkl)):
        P1 = scipy.linalg.expm(Q * float(childA.getBranchLength()))  
        P2 = scipy.linalg.expm(Q * float(childB.getBranchLength())) 
        vec1 = np.matmul(childA_lkl[pos_list], P1)
        vec2 = np.matmul(childB_lkl[pos_list], P2)
        vec_anc = [vec1[i] * vec2[i] for i in range(4)]
        ## !! CHECKPOINT !! ##
        #print("node",childA.parent.name,vec_anc) ## we are geiting 4 values/pos i.e lkl of all 4 nucleotides at each position.
        anc_AB_lkl.append(np.array(vec_anc)) ## position lkl to main list for this node

    return anc_AB_lkl

def calculate_lkl_recursively(node, Q, lkl_dict):
    if node.seq_leaf is not None:
        return # skip leaf nodes

    childA = node.children[0]
    childB = node.children[1]

    if childA.name not in lkl_dict.keys():
        calculate_lkl_recursively(childA, Q, lkl_dict)
    if childB.name not in lkl_dict.keys():
        calculate_lkl_recursively(childB, Q, lkl_dict)

    anc_name = node.name
    anc_AB_lkl = calculate_lkl_anc(childA, childB, Q, lkl_dict)
    lkl_dict[anc_name] = anc_AB_lkl

def main():
    parser = argparse.ArgumentParser(description="Process input files and calculate likelihoods.")
    parser.add_argument("--tree_path", type=str, default="data/tree.dat", help="Path to the tree topology file.")
    parser.add_argument("--branch_lengths_path", type=str, default="data/branchlength.dat", help="Path to the branch lengths file.")
    parser.add_argument("--msa_path", type=str, default="data/msa.dat", help="Path to the multiple sequence alignment file.")
    parser.add_argument("--mu", type=float, default=1, help="Value of mu (optional).")
    parser.add_argument("--output_lkl_root_node_only", action="store_true", help="Output likelihood of root node only.")

    args = parser.parse_args()

    tree_path = args.tree_path
    branch_lengths_path = args.branch_lengths_path
    msa_path = args.msa_path
    mu = args.mu
    output_lkl_root_node_only = args.output_lkl_root_node_only

    Q = np.array([[-3 * mu, mu, mu, mu], [mu, -3 * mu, mu, mu], [mu, mu, -3 * mu, mu], [mu, mu, mu, -3 * mu]])

    node_dict = read_input_files(tree_path, branch_lengths_path, msa_path)

    lkl_dict = {}

    for node in node_dict.values():
        if node.getSequence() is not None:
            lkl_dict[node.name] = node.get_lkl_list()

    parent_nodes_visited = []

    for root_node in node_dict.values():
        calculate_lkl_recursively(root_node, Q, lkl_dict)

    if output_lkl_root_node_only:
        npmm_sum=0
        root_name, root_value = list(lkl_dict.items())[-1]
        for pos_list in root_value:
            npmm=np.matmul(pos_list,[0.25,0.25,0.25,0.25])
            npmm_sum+=math.log(npmm)
        print(round(npmm_sum,4))
    else:
        for node in lkl_dict:
            npmm_sum=0
            for pos_list in lkl_dict[node]:
                npmm=np.matmul(pos_list,[0.25,0.25,0.25,0.25])
                npmm_sum+=math.log(npmm)
            print("Node ",node,"npmm: ",npmm_sum)

if __name__ == "__main__":
    main()
