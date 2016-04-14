#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *
# add - sourya
import Process_Queues
from Process_Queues import *

#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses one representative taxon of that taxa cluster
"""
def Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			"""
			clust_species_list[i] and clust_species_list[j]
			contain two taxa list of one or more elements
			"""
			entry = FindAvgXL(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 1)
			DistMat[j][i] = DistMat[i][j] = entry

	return
	
#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses aerage information of that taxa cluster
@param type_of_output: if 0, computes the average of XL measures
													1, returns the minimum of XL measures
													2, returns the maximum of XL measures
"""
def Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			# here both clust_species_list[i] and clust_species_list[j]
			# are one element lists (according to their construction)
			# we have extracted the corresponding element by using [0] operator (extracting first element)
			entry = FindAvgXL(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, 2)
			DistMat[j][i] = DistMat[i][j] = entry
	
	return

#-------------------------------------------
"""
following function checks if we are merging two leaves 
which do not have R3 as their consensus relation
in such a case, the function returns FALSE
"""
def CheckMergeImproperLeaves(clust_species_list, i, j):
	if (IsLeafCluster(clust_species_list, i) == True) and (IsLeafCluster(clust_species_list, j) == True):
		key1 = (clust_species_list[i][0], clust_species_list[j][0])
		key2 = (clust_species_list[j][0], clust_species_list[i][0])
		if key1 in TaxaPair_Reln_Dict:
			if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(RELATION_R3) == False):
				return True
		if key2 in TaxaPair_Reln_Dict:
			if (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(RELATION_R3) == False):
				return True
	
	return False

#-------------------------------------------
"""
this function finds a single minimum from the input matrix
"""
def Find_Unique_Min(Norm_DistMat, no_of_clust, clust_species_list):
	# initialize
	min_val = Norm_DistMat[0][1]
	min_idx_i = 0
	min_idx_j = 1
	# traverse through the matrix elements
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (i == j):
				continue
			
			if (Norm_DistMat[i][j] < min_val):
				if (CheckMergeImproperLeaves(clust_species_list, i, j) == False):
					min_val = Norm_DistMat[i][j]
					min_idx_i = i
					min_idx_j = j
			elif (Norm_DistMat[i][j] == min_val):
				if (CheckMergeImproperLeaves(clust_species_list, i, j) == False):
					if (GetR3RelnLevelConsVal(clust_species_list[i][0], clust_species_list[j][0]) > \
						GetR3RelnLevelConsVal(clust_species_list[min_idx_i][0], clust_species_list[min_idx_j][0])):
						min_idx_i = i
						min_idx_j = j
	
	return min_idx_i, min_idx_j

#--------------------------------------------------------
# this function is a shortcut to obtain the normalized expression 
# used in the agglomerative clustering proposed in this code
# as various methods are experimented, corresponding various forms of 
# agglomerative clustering is tried
#--------------------------------------------------------
def ObtainNormalizedVal(num, denom1, denom2):
  if ((denom1 + denom2) > 0):
    return (num * 1.0) / (denom1 + denom2)
  else:
    return 0

##---------------------------------------------
""" 
function to print the matrix content
N is the matrix dimension
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile):
  fp = open(textfile, 'a')
  fp.write('\n printing contents of ' + str(inp_str) + ' ---- ')
  for i in range(N):
    fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
    for j in range(i+1):
      fp.write(' ' + str(inp_data[i][j]))
  fp.close()

#-----------------------------------------------------------------
def GetR3RelnLevelConsVal(x1, x2):
	key1 = (x1, x2)
	key2 = (x2, x1)
	if key1 in TaxaPair_Reln_Dict:	
		return TaxaPair_Reln_Dict[key1]._GetRelnLevelDiff(RELATION_R3, -1, False, True)
	elif key2 in TaxaPair_Reln_Dict:
		return TaxaPair_Reln_Dict[key2]._GetRelnLevelDiff(RELATION_R3, -1, False, True)

#-------------------------------------------
"""
this function processes input distance matrix in every iteration
and finds the pair of indices satisfying minimum distance criterion 
used in NJ based algorithm
"""
def Get_NJ_Based_Min_Pair_Idx(DistMat, Norm_DistMat, no_of_clust, clust_species_list, \
	NJ_RULE_USED, Output_Text_File):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n\n\n\n **** iteration start --- number of clusters: ' + str(no_of_clust))
		fp.write('\n clust_species_list : ' + str(clust_species_list))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, DistMat, 'DistMat', Output_Text_File)
	
	# allocate one new square matrix which will contain the 
	# normalized matrix elements (w.r.t the sum of sum of rows and columns)
	sum_list = []
	for i in range(no_of_clust):
		t = 0
		for j in range(no_of_clust):
			t = t + DistMat[i][j]
		sum_list.append(t)

	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (NJ_RULE_USED == AGGLO_CLUST):
				# we normalize the extra lineage based score
				# by the sum of extra lineages for all other taxa from the taxa indices i and j
				# modified - sourya
				#Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
				# add - sourya
				Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i] - DistMat[i][j], \
					sum_list[j] - DistMat[i][j])
				# modified - sourya
				#Norm_DistMat[i][j] = DistMat[i][j]
				# end add - sourya
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
			else:
				ri = sum_list[i] / (no_of_clust - 2)
				rj = sum_list[j] / (no_of_clust - 2)
				Norm_DistMat[i][j] = (DistMat[i][j] - ri - rj)
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
				
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n printing contents of sum_list --- ' + str(sum_list))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, Norm_DistMat, 'Norm_DistMat', Output_Text_File)

	# now we have to find the minimum among these elements 
	# present in the matrix Norm_DistMat
	min_idx_i, min_idx_j = Find_Unique_Min(Norm_DistMat, no_of_clust, clust_species_list)

	return min_idx_i, min_idx_j

#-------------------------------------------
"""
checks whether a taxa cluster specified by the input index is a leaf
"""
def IsLeafCluster(clust_species_list, idx):
	if (len(clust_species_list[idx]) == 1):
		return True
	return False

#-------------------------------------------
"""
this function places one subtree at an edge of a second subtree
also adjusts their parent information
parameters:
1) source subtree (which needs to be re-positioned)
2) edge of destination subtree - indicated by the child node
3) Curr_tree: Tree containing all these subtrees

It creates one new internal node within the edge of the destination subtree
and places the source subtree as its children
"""
def InsertSubTree(Curr_tree, src_subtree, child_dest_subtree, Output_Text_File):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of src_subtree: ' + str(Node_Label(src_subtree)))      
		fp.write('\n label of child_dest_subtree: ' + str(Node_Label(child_dest_subtree)))
		fp.close()
	
	# create new internal node 
	newnode = Node()
	
	parent_dest_subtree = child_dest_subtree.parent_node
	parent_src_subtree = src_subtree.parent_node
	
	# its parent node will be the previous MRCA node of all the taxa in two clusters
	parent_dest_subtree.add_child(newnode)
	newnode.parent_node = parent_dest_subtree
	parent_dest_subtree.remove_child(child_dest_subtree)
	child_dest_subtree.parent_node = None
	newnode.add_child(child_dest_subtree)
	child_dest_subtree.parent_node = newnode
	
	if (parent_src_subtree is not None):
		parent_src_subtree.remove_child(src_subtree)
		src_subtree.parent_node = None
	newnode.add_child(src_subtree)
	src_subtree.parent_node = newnode
	
	# update splits of the resulting tree
	Curr_tree.update_splits(delete_outdegree_one=False)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of dest_subtree: ' + str(Node_Label(child_dest_subtree.parent_node)))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function has following parameters:
1) first_cluster_mrca_node: root of 1st subtree 
2) second_cluster_mrca_node: root of 2nd subtree 
3) all_taxa_mrca_node: root of all these trees
4) Curr_tree: Tree containing all these subtrees

It creates one new internal node as a child of all_taxa_mrca_node
and places above mentioned subtrees as its children
"""
def MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, \
	all_taxa_mrca_node, taxa_list, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
		fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
		fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
		fp.close()

	# create new internal node 
	newnode = Node()
	
	# its parent node will be the previous MRCA node of all the taxa in two clusters
	all_taxa_mrca_node.add_child(newnode)
	newnode.parent_node = all_taxa_mrca_node
	all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
	first_cluster_mrca_node.parent_node = None
	all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
	second_cluster_mrca_node.parent_node = None
	
	# add these individual clusters' MRCA node as its children
	newnode.add_child(first_cluster_mrca_node)
	first_cluster_mrca_node.parent_node = newnode
	newnode.add_child(second_cluster_mrca_node)
	second_cluster_mrca_node.parent_node = newnode
	
	# update splits of the resulting tree
	Curr_tree.update_splits(delete_outdegree_one=False)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of all taxa mrca node (recomputed): ' \
			+ str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function merges two leaf nodes in the non-refined species tree
"""
def Merge_Leaves(Curr_tree, clust_species_list, idx1, idx2, taxa_list, Output_Text_File):
	"""
	both clusters are leaves - check whther there exist any level difference of them
	otherwise they can be treated as a sibling taxa pair in the species tree
	"""
	taxa1 = clust_species_list[idx1][0]
	taxa2 = clust_species_list[idx2][0]

	clust1_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx1])
	clust2_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx2])
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, \
		all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree


# ---------------------------------------
"""
this function checks whether taxa1 can be related with taxa2 via the target_reln
given their relation frequencies and also the level values
"""
def Freq_Level_Higher(taxa1, taxa2, target_reln):
	key1 = (taxa1, taxa2)
	key2 = (taxa2, taxa1)
	if key1 in TaxaPair_Reln_Dict:
		return TaxaPair_Reln_Dict[key1].CheckHigherPriority(target_reln)
	elif key2 in TaxaPair_Reln_Dict:
		return TaxaPair_Reln_Dict[key2].CheckHigherPriority(Complementary_Reln(target_reln))

	return 0

#-------------------------------------------
"""
this function merges one leaf node and another non leaf node (taxa cluster)
to refine the output species tree
"""
def Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, leaf_idx, non_leaf_idx, \
	taxa_list, Output_Text_File, DIST_MAT_TYPE):
	
	#----------------------------------------------------------------
	# representing first cluster as A (leaf node) and second cluster as (B,C)
	#----------------------------------------------------------------
	
	# this is the MRCA node corresponding to the input leaf node
	leaf_taxon = clust_species_list[leaf_idx][0]
	leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[leaf_idx])
	
	# this is the MRCA node corresponding to the non leaf taxa cluster
	non_leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[non_leaf_idx])
	
	# MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	clust2_children = non_leaf_mrca_node.child_nodes()
	clust2_child1_taxa_list = []
	clust2_child2_taxa_list = []
	
	clust2_child1_taxa_list = GetPreorderTaxaList(clust2_children[0], clust2_child1_taxa_list, taxa_list)
	clust2_child2_taxa_list = GetPreorderTaxaList(clust2_children[1], clust2_child2_taxa_list, taxa_list)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging one leaf and other non leaf node') 
		fp.write('\n ---- First leaf idx list: ' + str(clust_species_list[leaf_idx]))
		fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[non_leaf_idx]))
		fp.write('\n ---- clust2_child1_taxa_list: ' + str(clust2_child1_taxa_list)) 
		fp.write('\n ---- clust2_child2_taxa_list: ' + str(clust2_child2_taxa_list))
		fp.close()
	
	#----------------------------------------------------

	# obtain three sets of XL values
	"""
	depending on the cardinality of constituent taxa list
	we employ cluster based average or single taxa pair based XL
	"""
	if (len(clust2_child1_taxa_list) == 1) and (len(clust2_child2_taxa_list) == 1):
		XL_clust2_child1_clust2_child2 = FindAvgXL(clust2_child1_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust2_child1_clust2_child2 = FindAvgXL(clust2_child1_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 1)
		
	if (len(clust2_child1_taxa_list) == 1):
		XL_clust2_child1_leaf = FindAvgXL(clust2_child1_taxa_list, clust_species_list[leaf_idx], \
			DIST_MAT_TYPE, 0)
	else:
		XL_clust2_child1_leaf = FindAvgXL(clust2_child1_taxa_list, clust_species_list[leaf_idx], \
			DIST_MAT_TYPE, 1)
		
	if (len(clust2_child2_taxa_list) == 1):
		XL_clust2_child2_leaf = FindAvgXL(clust2_child2_taxa_list, clust_species_list[leaf_idx], \
			DIST_MAT_TYPE, 0)
	else:
		XL_clust2_child2_leaf = FindAvgXL(clust2_child2_taxa_list, clust_species_list[leaf_idx], \
			DIST_MAT_TYPE, 1)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ---- XL_clust2_child1_clust2_child2: ' + str(XL_clust2_child1_clust2_child2))
		fp.write('\n ---- XL_clust2_child1_leaf: ' + str(XL_clust2_child1_leaf))
		fp.write('\n ---- XL_clust2_child2_leaf: ' + str(XL_clust2_child2_leaf))
		fp.close()

	XL_list = [XL_clust2_child1_clust2_child2, XL_clust2_child1_leaf, XL_clust2_child2_leaf]
	
	#---------------------------------------------
	"""
	Note: A is the external leaf
	(B,C) is the existing cluster
	
	default configuration is (A,(B,C)) 
	however, we explore other options as well
	"""
	
	# construct single element lists from each representative class
	A_Single_Taxa_List = []
	A_Single_Taxa_List.append(leaf_taxon)
	B_Single_Taxa_List = []
	B_Single_Taxa_List.append(clust2_child1_taxa_list[0])
	C_Single_Taxa_List = []
	C_Single_Taxa_List.append(clust2_child2_taxa_list[0])
	
	"""
	case 1.1 - check if B can be placed externally, and (A,C) together can be placed as siblings
	this can be done only if C is not a leaf
	"""
	# modified - sourya
	if (len(clust2_child2_taxa_list) > 1):
		if ((len(clust2_child1_taxa_list) > 1) and (CheckR1Reln(B_Single_Taxa_List, A_Single_Taxa_List) > 0)) \
			or ((len(clust2_child1_taxa_list) == 1) \
				and (Freq_Level_Higher(B_Single_Taxa_List[0], A_Single_Taxa_List[0], RELATION_R1) == 1)):
				
	#if (len(clust2_child2_taxa_list) > 1) and (CheckR1Reln(B_Single_Taxa_List, A_Single_Taxa_List) > 0):	# comment - sourya
			"""
			if R1(B,A) will be predominant compared to R1(A,B)
			then the configuration (B, (A,C)) can be employed
			"""
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 1 --- Merging Condition (B, (A, C))')
				fp.close()
			src_subtree_node = leaf_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree

	"""
	case 1.2 - check if C can be placed externally, and (A,B) together can be placed as siblings
	this can be done only if B is not a leaf
	"""
	# modified - sourya
	if (len(clust2_child1_taxa_list) > 1):
		if ((len(clust2_child2_taxa_list) > 1) and (CheckR1Reln(C_Single_Taxa_List, A_Single_Taxa_List) > 0)) \
			or ((len(clust2_child2_taxa_list) == 1) \
				and (Freq_Level_Higher(C_Single_Taxa_List[0], A_Single_Taxa_List[0], RELATION_R1) == 1)):
				
	#if (len(clust2_child1_taxa_list) > 1) and (CheckR1Reln(C_Single_Taxa_List, A_Single_Taxa_List) > 0):	# comment - sourya
			"""
			if R1(C,A) will be predominant compared to R1(A,C)
			then the configuration (C, (A,B)) can be employed
			"""
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 2 --- Merging Condition (C, (A, B))')
				fp.close()
			src_subtree_node = leaf_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree

	#---------------------------------------------
	# XL based check

	if (XL_clust2_child1_clust2_child2 == min(XL_list)):
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n *** Case 3 --- Merging Condition (A, (B, C))')
			fp.close()
		# configuration (A,(B,C)) will be employed
		Curr_tree = MergeSubtrees(Curr_tree, leaf_mrca_node, non_leaf_mrca_node, \
			all_taxa_mrca_node, taxa_list, Output_Text_File)
		return Curr_tree
	elif (XL_clust2_child2_leaf == min(XL_list)):
		if (len(clust2_child2_taxa_list) > 1):
			# configuration (B, (A, C)) will be employed
			# provided C is not leaf
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 4 --- Merging Condition (B, (A, C))')
				fp.close()
			src_subtree_node = leaf_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
	elif (XL_clust2_child1_leaf == min(XL_list)):
		if (len(clust2_child1_taxa_list) > 1):
			# configuration (C, (A, B)) will be employed 
			# provided B is not leaf
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 5 --- Merging Condition (C, (A, B))')
				fp.close()
			src_subtree_node = leaf_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
	
	# otherwise configuration (A,(B,C)) will be employed
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n *** Case 6 --- Merging Condition (A, (B, C))')
		fp.close()
	Curr_tree = MergeSubtrees(Curr_tree, leaf_mrca_node, non_leaf_mrca_node, \
		all_taxa_mrca_node, taxa_list, Output_Text_File)
	return Curr_tree

#-------------------------------------------
"""
this function merges both non leaf nodes (taxa cluster) to refine the output species tree
"""
def Merge_Both_NonLeaf(Curr_tree, clust_species_list, idx1, idx2, taxa_list, \
	Output_Text_File, DIST_MAT_TYPE):
	
	#----------------------------------------------------------------
	#representing first cluster as (A,B) and second cluster as (C,D)
	#----------------------------------------------------------------
	
	# this is the MRCA node corresponding to the first cluster taxa list
	clust1_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx1])
	
	# this is the MRCA node corresponding to the second cluster taxa list
	clust2_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx2])
	
	# MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	# child nodes of both these MRCA nodes of respective taxa clusters
	clust1_children = clust1_mrca_node.child_nodes()
	clust1_child1_taxa_list = []
	clust1_child2_taxa_list = []

	clust1_child1_taxa_list = GetPreorderTaxaList(clust1_children[0], clust1_child1_taxa_list, taxa_list)
	clust1_child2_taxa_list = GetPreorderTaxaList(clust1_children[1], clust1_child2_taxa_list, taxa_list)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging both non leaf nodes') 
		fp.write('\n ---- First non leaf idx list: ' + str(clust_species_list[idx1]))
		fp.write('\n ---- clust1_child1_taxa_list: ' + str(clust1_child1_taxa_list)) 
		fp.write('\n ---- clust1_child2_taxa_list: ' + str(clust1_child2_taxa_list))
		fp.close()
	
	clust2_children = clust2_mrca_node.child_nodes()
	clust2_child1_taxa_list = []
	clust2_child2_taxa_list = []

	clust2_child1_taxa_list = GetPreorderTaxaList(clust2_children[0], clust2_child1_taxa_list, taxa_list)
	clust2_child2_taxa_list = GetPreorderTaxaList(clust2_children[1], clust2_child2_taxa_list, taxa_list)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging both non leaf nodes') 
		fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[idx2]))
		fp.write('\n ---- clust2_child1_taxa_list: ' + str(clust2_child1_taxa_list)) 
		fp.write('\n ---- clust2_child2_taxa_list: ' + str(clust2_child2_taxa_list))
		fp.close()
	
	#---------------------------------------------------------
	# obtain XL values for different taxa sets
	# XL(A,B)
	if (len(clust1_child1_taxa_list) == 1) and (len(clust1_child2_taxa_list) == 1):
		XL_clust1_child1_clust1_child2 = FindAvgXL(clust1_child1_taxa_list, \
			clust1_child2_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust1_child1_clust1_child2 = FindAvgXL(clust1_child1_taxa_list, \
			clust1_child2_taxa_list, DIST_MAT_TYPE, 1)
		
	# XL(C,D)
	if (len(clust2_child1_taxa_list) == 1) and (len(clust2_child2_taxa_list) == 1):
		XL_clust2_child1_clust2_child2 = FindAvgXL(clust2_child1_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust2_child1_clust2_child2 = FindAvgXL(clust2_child1_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 1)
		
	# XL(A,C)
	if (len(clust1_child1_taxa_list) == 1) and (len(clust2_child1_taxa_list) == 1):
		XL_clust1_child1_clust2_child1 = FindAvgXL(clust1_child1_taxa_list, \
			clust2_child1_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust1_child1_clust2_child1 = FindAvgXL(clust1_child1_taxa_list, \
			clust2_child1_taxa_list, DIST_MAT_TYPE, 1)
	
	# XL(A,D)
	if (len(clust1_child1_taxa_list) == 1) and (len(clust2_child2_taxa_list) == 1):
		XL_clust1_child1_clust2_child2 = FindAvgXL(clust1_child1_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust1_child1_clust2_child2 = FindAvgXL(clust1_child1_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 1)
	
	# XL(B,C)
	if (len(clust1_child2_taxa_list) == 1) and (len(clust2_child1_taxa_list) == 1):
		XL_clust1_child2_clust2_child1 = FindAvgXL(clust1_child2_taxa_list, \
			clust2_child1_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust1_child2_clust2_child1 = FindAvgXL(clust1_child2_taxa_list, \
			clust2_child1_taxa_list, DIST_MAT_TYPE, 1)
	
	# XL(B,D)
	if (len(clust1_child2_taxa_list) == 1) and (len(clust2_child2_taxa_list) == 1):
		XL_clust1_child2_clust2_child2 = FindAvgXL(clust1_child2_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 0)
	else:
		XL_clust1_child2_clust2_child2 = FindAvgXL(clust1_child2_taxa_list, \
			clust2_child2_taxa_list, DIST_MAT_TYPE, 1)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ---- XL_clust1_child1_clust1_child2: ' + str(XL_clust1_child1_clust1_child2))
		fp.write('\n ---- XL_clust2_child1_clust2_child2: ' + str(XL_clust2_child1_clust2_child2))
		fp.write('\n ---- XL_clust1_child1_clust2_child1: ' + str(XL_clust1_child1_clust2_child1))
		fp.write('\n ---- XL_clust1_child1_clust2_child2: ' + str(XL_clust1_child1_clust2_child2))
		fp.write('\n ---- XL_clust1_child2_clust2_child1: ' + str(XL_clust1_child2_clust2_child1))
		fp.write('\n ---- XL_clust1_child2_clust2_child2: ' + str(XL_clust1_child2_clust2_child2))
		fp.close()
	
	if (XL_clust1_child1_clust1_child2 <= XL_clust2_child1_clust2_child2):
		
		"""
		cluster (A,B) will remain intact
		we have to inspect between three configurations
		1) ((A,B),(C,D)), 2) (C, (D, (A, B))), 3) (D, (C, (A, B)))
		"""
		
		#-------------------------------------
		# add - sourya
		# condition if both C and D are taxa clusters (length > 1)
		if (len(clust2_child1_taxa_list) > 1) and (len(clust2_child2_taxa_list) > 1):
			C_Single_Taxa_List = []
			C_Single_Taxa_List.append(clust2_child1_taxa_list[0])
			D_Single_Taxa_List = []
			D_Single_Taxa_List.append(clust2_child2_taxa_list[0])
			AB_Combined_Taxa_List = []
			AB_Combined_Taxa_List.append(clust_species_list[idx1][0])
			"""
			we check if both C and D are at least levels higher than AB
			"""
			if (CheckR1Reln(C_Single_Taxa_List, AB_Combined_Taxa_List, False) > 0) \
				and (CheckR1Reln(D_Single_Taxa_List, AB_Combined_Taxa_List, False) > 0):
				""" 
				we check the XL(C,AB) and XL(D,AB)
				whichever is lower, assign AB as a sibling to that cluster
				"""
				XL_C_AB = max(XL_clust1_child1_clust2_child1, XL_clust1_child2_clust2_child1)
				XL_D_AB = max(XL_clust1_child1_clust2_child2, XL_clust1_child2_clust2_child2)
				if (XL_C_AB < XL_D_AB):
					"""
					configuration (D, (C, (A, B))) can be employed here
					"""
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n *** Case 1A - Merging Condition (D, (C, (A, B)))')
						fp.close()
					src_subtree_node = clust1_mrca_node
					dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
					Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
					return Curr_tree
				else:
					"""
					configuration (C, (D, (A, B))) can be employed here
					"""
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n *** Case 1B - Merging Condition (C, (D, (A, B)))')
						fp.close()
					src_subtree_node = clust1_mrca_node
					dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
					Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
					return Curr_tree
		
		# end add - sourya
		#-------------------------------------
		
		#  XL(C,:) / 2
		# modify - sourya
		# replacing average with maximum operator
		#XL_clust2_child1_clust1_avg = (XL_clust1_child1_clust2_child1 + XL_clust1_child2_clust2_child1) / 2.0
		XL_clust2_child1_clust1_avg = max(XL_clust1_child1_clust2_child1, XL_clust1_child2_clust2_child1)
		# XL(D,:) / 2
		# modify - sourya
		# replacing average with maximum operator
		#XL_clust2_child2_clust1_avg = (XL_clust1_child1_clust2_child2 + XL_clust1_child2_clust2_child2) / 2.0
		XL_clust2_child2_clust1_avg = max(XL_clust1_child1_clust2_child2, XL_clust1_child2_clust2_child2)
		
		"""
		# case 1A - XL(C,D) is lower than both XL(C,:) / 2 and XL(D,:) / 2 (: denotes A and B)
		"""
		if (XL_clust2_child1_clust2_child2 <= XL_clust2_child1_clust1_avg) and (XL_clust2_child1_clust2_child2 <= XL_clust2_child2_clust1_avg):
			# ((A,B),(C,D)) configuration can be employed
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 3 - Default Merging Condition ((A,B),(C,D))')
				fp.close()
			Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
			return Curr_tree
		
		"""
		# case 2A - XL(C,:) / 2 is lower than both XL(C,D) and XL(D,:) / 2 (: denotes A and B)
		"""
		if (XL_clust2_child1_clust1_avg <= XL_clust2_child1_clust2_child2) and (XL_clust2_child1_clust1_avg <= XL_clust2_child2_clust1_avg):
			# configuration (D, (C, (A, B))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 4 - Merging Condition (D, (C, (A, B)))')
				fp.close()
			src_subtree_node = clust1_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
		"""
		# case 3A - XL(D,:) / 2 is lower than both XL(C,D) and XL(C,:) / 2 (: denotes A and B)
		"""
		if (XL_clust2_child2_clust1_avg <= XL_clust2_child1_clust1_avg) and (XL_clust2_child2_clust1_avg <= XL_clust2_child1_clust2_child2):
			# configuration (C, (D, (A, B))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 5 - Merging Condition (C, (D, (A, B)))')
				fp.close()
			src_subtree_node = clust1_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
	else:
		
		"""
		cluster (C,D) will remain intact
		we have to inspect between three configurations
		1) ((A,B),(C,D)), 2) (A, (B, (C, D))), 3) (B, (A, (C, D)))
		"""

		#-------------------------------------
		# add - sourya
		# condition if both A and B are taxa clusters (length > 1)
		if (len(clust1_child1_taxa_list) > 1) and (len(clust1_child2_taxa_list) > 1):
			A_Single_Taxa_List = []
			A_Single_Taxa_List.append(clust1_child1_taxa_list[0])
			B_Single_Taxa_List = []
			B_Single_Taxa_List.append(clust1_child2_taxa_list[0])
			CD_Combined_Taxa_List = []
			CD_Combined_Taxa_List.append(clust_species_list[idx2][0])
			"""
			we check if both A and B are at least levels higher than CD
			"""
			if (CheckR1Reln(A_Single_Taxa_List, CD_Combined_Taxa_List, False) > 0) \
				and (CheckR1Reln(B_Single_Taxa_List, CD_Combined_Taxa_List, False) > 0):
				""" 
				we check the XL(A,CD) and XL(B,CD)
				whichever is lower, assign CD as a sibling to that cluster
				"""
				XL_A_CD = max(XL_clust1_child1_clust2_child1, XL_clust1_child1_clust2_child2)
				XL_B_CD = max(XL_clust1_child2_clust2_child1, XL_clust1_child2_clust2_child2)
				if (XL_A_CD < XL_B_CD):
					"""
					the configuration (B, (A, (C, D))) can be employed here
					"""
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n *** Case 2A - Merging Condition (B, (A, (C, D)))')
						fp.close()
					src_subtree_node = clust2_mrca_node
					dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child1_taxa_list)
					Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
					return Curr_tree
				else:
					"""
					configuration (A, (B, (C, D))) can be employed here
					"""
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n *** Case 10 - Merging Condition (A, (B, (C, D)))')
						fp.close()
					src_subtree_node = clust2_mrca_node
					dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child2_taxa_list)
					Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
					return Curr_tree
		
		# end add - sourya
		#-------------------------------------

		#  XL(A,:) / 2
		# modify - sourya
		# replacing average with maximum operator
		#XL_clust1_child1_clust2_avg = (XL_clust1_child1_clust2_child1 + XL_clust1_child1_clust2_child2) / 2.0
		XL_clust1_child1_clust2_avg = max(XL_clust1_child1_clust2_child1, XL_clust1_child1_clust2_child2)
		# XL(B,:) / 2
		# modify - sourya
		# replacing average with maximum operator
		#XL_clust1_child2_clust2_avg = (XL_clust1_child2_clust2_child1 + XL_clust1_child2_clust2_child2) / 2.0
		XL_clust1_child2_clust2_avg = max(XL_clust1_child2_clust2_child1, XL_clust1_child2_clust2_child2)
		"""
		# case 1B - XL(A,B) is lower than both XL(A,:) / 2 and XL(B,:) / 2 (: denotes C and D)
		"""
		if (XL_clust1_child1_clust1_child2 <= XL_clust1_child1_clust2_avg) \
			and (XL_clust1_child1_clust1_child2 <= XL_clust1_child2_clust2_avg):
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 8 - Default Merging Condition ((A,B),(C,D))')
				fp.close()
			Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, \
				all_taxa_mrca_node, taxa_list, Output_Text_File)
			return Curr_tree
		
		"""
		# case 2B - XL(A,:) / 2 is lower than both XL(A,B) and XL(B,:) / 2 (: denotes C and D)
		"""
		if (XL_clust1_child1_clust2_avg <= XL_clust1_child1_clust1_child2) \
			and (XL_clust1_child1_clust2_avg <= XL_clust1_child2_clust2_avg):
			# the configuration (B, (A, (C, D))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 9 - Merging Condition (B, (A, (C, D)))')
				fp.close()
			src_subtree_node = clust2_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
			
		"""
		# case 3B - XL(B,:) / 2 is lower than both XL(A,B) and XL(A,:) / 2 (: denotes C and D)
		"""
		if (XL_clust1_child2_clust2_avg <= XL_clust1_child1_clust1_child2) \
			and (XL_clust1_child2_clust2_avg <= XL_clust1_child1_clust2_avg):
			# the configuration (A, (B, (C, D))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Case 10 - Merging Condition (A, (B, (C, D)))')
				fp.close()
			src_subtree_node = clust2_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
			
	# default merging if none of the above conditions work ((A,B),(C,D)) configuration can be employed
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n *** Case 11 - Default Merging Condition ((A,B),(C,D))')
		fp.close()
	Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, \
		all_taxa_mrca_node, taxa_list, Output_Text_File)
	return Curr_tree

#-------------------------------------------
"""
this is a new function to merge taxa clusters for agglomerative clustering - sourya
"""
def Merge_Cluster_Pair_New(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, \
	Output_Text_File, DIST_MAT_TYPE):

	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	if (isleaf_clust1 == True) and (isleaf_clust2 == True):
		Curr_tree = Merge_Leaves(Curr_tree, clust_species_list, min_idx_i, min_idx_j, \
			taxa_list, Output_Text_File)
	elif (isleaf_clust1 == True) and (isleaf_clust2 == False):
		Curr_tree = Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, \
			Output_Text_File, DIST_MAT_TYPE)
	elif (isleaf_clust1 == False) and (isleaf_clust2 == True):
		Curr_tree = Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, min_idx_j, min_idx_i, taxa_list, \
			Output_Text_File, DIST_MAT_TYPE)
	elif (isleaf_clust1 == False) and (isleaf_clust2 == False):
		Curr_tree = Merge_Both_NonLeaf(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, \
			Output_Text_File, DIST_MAT_TYPE)
	
	return Curr_tree

#-------------------------------------------
"""
this function merges a pair of clusters whose indices are pointed by the min_idx_i and min_idx_j entries
this is part of the proposed agglomerative clustering
taxa_list is the union of these two clusters (species contents)
"""
def Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File):
	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	
	if (isleaf_clust1):
		first_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
	else:
		first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
	
	if (isleaf_clust2):
		second_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
	else:
		second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
	
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	Curr_tree = MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, \
		all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree

#-------------------------------------------
"""
this function processes one internal node (basically the children list)
to resolve multifurcation
"""
def ResolveMultifurcation(Curr_tree, clust_species_list, no_of_input_clusters, Output_Text_File, \
	NJ_RULE_USED, DIST_MAT_TYPE, NJ_MERGE_CLUST):
	# total number of clusters
	no_of_clust = no_of_input_clusters

	# allocate a 2D square matrix of no_of_clust dimension
	DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
	Norm_DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)

	# here we compute the ILS score of the current cluster pair
	# with respect to input gene tree list
	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n Examining ILS score for individual cluster pairs ')
		fp.close()      

	#---------------------------------------
	## using single taxon as a representative of the taxa cluster
	#Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE)

	# using average information from a taxa cluster
	Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE)
	#---------------------------------------

	# loop to execute the agglomerative clustering
	while(no_of_clust > 2):
		min_idx_i, min_idx_j = Get_NJ_Based_Min_Pair_Idx(DistMat, Norm_DistMat, \
			no_of_clust, clust_species_list, NJ_RULE_USED, Output_Text_File)

		# note down the taxa list in these two indices of the clust_species_list
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete species list ' + str(taxa_list))
			fp.close()
			
		"""
		now we merge the pair of clusters pointed by these indices
		"""
		if (NJ_MERGE_CLUST == 1):
			Curr_tree = Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, \
				taxa_list, Output_Text_File)
		else:
			Curr_tree = Merge_Cluster_Pair_New(Curr_tree, clust_species_list, \
				min_idx_i, min_idx_j, taxa_list, Output_Text_File, DIST_MAT_TYPE)
		#---------------------------------------------------------      
		# remove individual clusters' taxa information from the clust_species_list
		# and add taxa_list as a new element
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		
		# comment - sourya
		#clust_species_list.append(taxa_list)    
		
		# add - sourya
		taxa_list_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
		preorder_taxa_list = []
		for n in taxa_list_mrca_node.preorder_iter():
			if (n.is_leaf() == True):
				if n.taxon.label in taxa_list:
					preorder_taxa_list.append(n.taxon.label)
		clust_species_list.append(preorder_taxa_list)   
		# end add - sourya
		#---------------------------------------------------------      
		"""
		adjust the DistMat by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		# first append one row
		DistMat = numpy.vstack((DistMat, numpy.zeros((1, no_of_clust), dtype=numpy.float)))
		# then append one column
		DistMat = numpy.hstack((DistMat, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.float)))
		# now apply reshape operation to get proper square matrix dimension
		DistMat = numpy.reshape(DistMat, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
		
		# now fill the elements of the new added row and column
		for k in range(no_of_clust):
			if (NJ_RULE_USED == AGGLO_CLUST):
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j]) / 2.0
			else:
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j] \
					- DistMat[min_idx_i][min_idx_j]) / 2.0
			# symmetric property
			DistMat[no_of_clust][k] = DistMat[k][no_of_clust]
		
		# now remove the rows and columns corresponding to min_idx_i and min_idx_j
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=1)	# delete the column
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=1)	# delete the column

		# clear Norm_DistMat
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=0)	# delete the row
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=1)	# delete the column
		Norm_DistMat.fill(0)
		
		# decrement the number of clusters considered
		no_of_clust = no_of_clust - 1
	
	# add - sourya
	"""
	delete the distance matrices
	"""
	del DistMat
	del Norm_DistMat
	# end add - sourya
	
	return
            
#-------------------------------------------
# this function refines input supertree such that the supertree becomes binary
# this is required for proper benchmarking with existing binary tree construction methods on 
# ILS sorting
def Refine_Supertree_Binary_Form(Curr_tree, Output_Text_File, NJ_RULE_USED, \
	DIST_MAT_TYPE, NJ_MERGE_CLUST):
	"""
	we traverse input tree internal nodes in postorder fashion
	and list the child nodes of it
	if the no of children > 2 then it is a case of multifurcation
	for resolving
	"""
	for curr_node in Curr_tree.postorder_internal_node_iter():
		curr_node_children = curr_node.child_nodes()
		if (len(curr_node_children) > 2):
			# create a list which will contain the species list lying under 
			# individual child nodes of rhe current node
			clust_species_list = []
			
			for x in curr_node_children:
				subl = []
				for n in x.preorder_iter():
					if (n.is_leaf() == True):
						subl.append(n.taxon.label)
				clust_species_list.append(subl)
			
			# call the resolving routine
			ResolveMultifurcation(Curr_tree, clust_species_list, len(curr_node_children), Output_Text_File, \
				NJ_RULE_USED, DIST_MAT_TYPE, NJ_MERGE_CLUST)

	return
