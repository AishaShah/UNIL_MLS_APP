### PSEUDOCODE ###

'''
Input: 
    1. tree in table format
    2. Branch lengths
    3. MSA (Multiple Sequence Alignment)
    4. Model of evolution <-- optional 
Output: 
    Probability of alignment
Creating a class named tree_class to store tree
    Each Node should be an object in tree_class
        Each object should have following attributes
            1. Branch lengths
            2. Sequence  of leaf nodes
            3. pointers to childs and parent i.e other objects

    Matrix <- pairwise alignment of all objects/sequences <-- instead of dict??

    Class Functions:
        get_branch_length
        get_sequence
        get_children
        get_nuc_list(): GET LIST OF NUCLEOTIDE PRESENT AT EACH SITE IN FORMAT [1,0,0,0] 
'''
#MAIN PROGRAMME
#we should have branch lengths, sequences (not all nodes ahve sequences--> only the leaf nodes) , tree in our class data structure

# this dict will store each node and the probabilities of 4 nucleotides at each position
# dict = {key is node name: values are list of lists}
lkl_dict =  {node_name:nuc_likelihood}   

for node in tree.get_all_nodes:
    if node.get_sequence() is not null: # if it is not internal node 
        nuc_list=node.get_nuc_list() # which nuc is present at each pposition --> a list of lists like [[1,0,0,0],[0,1,0,0]....]
        Probability_Matrix = scipy.linalg.expm(Q * node.get_branch_length()) # same for all positions
        for pos_list in nuc_list: # for each pos we have a list like [1,0,0,0] in nuc_list
            nuc_lkl = numpy.matmul(Probability_Matrix, pos_list ) # calculate liklihood of each nucleotide
            merge nuc_lkl to lkl_dict[node] # add to dict
   
    # if it is an internal node:
    else:
        child_names=node.get_children() # get_children will return a list of child node names --> two children per parent
        childA_lkl = lkl_dict[child_names[0]] # we will get a list of lists with prob for each nuclleotide at each position
        childB_lkl = lkl_dict[child_names[1]]
        for pos_list  in range(len(childA_lkl)): #using all positions in both children nodes, calculate ancetral nucleotide liklihoods 
            joint_lkl_AB = [childA_lkl[pos_list ][i] * childB_lkl[pos_list ][i] for i in range(4)]
            anc_lkl= numpy.matmul(joint_lkl_AB, [0.25,0.25,0.25,0.25] ) # all nucleotides ahve sam eprobabilities
            merge anc_lkl to lkl_dict[node]
return sum(math.log(anc_lkl)) # values in npmm will be the ones for root

