#!/usr/bin/env python

import Header
from Header import *

#-------------------------------------------
#"""
#this function prints the tree in Newick format
#this function is used in COSPEDSpec
#"""
#def PrintNewick(root_clust_node_idx):
	#if 0:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList()

	#Tree_Str_List = ''
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		#"""
		#set the explored status of the current node to true
		#"""
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		#"""
		#this is the list of taxa of this cluster
		#"""
		#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		#"""
		#get the out edge list of the current cluster which are not explored yet 
		#"""
		#outnodes = []
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#outnodes.append(l)
		
		#""" 
		#at first, print the contents of this taxa cluster
		#if the cluster has more than one taxon, then use ( and ) to enclose the taxa list
		#"""
		#if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + '('
			
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + '('
		#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
		#if (len(spec_list) > 1):
			#Tree_Str_List = Tree_Str_List + ')'
		#"""
		#here we check if the cluster has one or more out edges
		#then recursively traverse all the out edge clusters
		#"""
		#if (len(outnodes) > 0):
			## first add one comma
			#Tree_Str_List = Tree_Str_List + ','
			
			## then add one opening bracket, within which, all the out edge cluster contents will reside
			#Tree_Str_List = Tree_Str_List + '('
			
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					#if (i < (len(outnodes) - 1)):
						## we check whether any subsequent node belonging to the outnodes list
						## is left for traverse
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						## in this case, we append one comma
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','      
			
			## at last, append one closing bracket, signifying the end of out edge cluster contents
			#Tree_Str_List = Tree_Str_List + ')'

		#if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			#Tree_Str_List = Tree_Str_List + ')'
		
	#return Tree_Str_List    

#--------------------------------------------------------
"""
this is a modified function - sourya - April 7, 2016
this function prints the tree in Newick format
"""
def PrintNewick(root_clust_node_idx):
	if 0:
		print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList()

	"""
	initialize the tree representative newick string
	"""
	Tree_Str_List = ''
	
	"""
	process the node provided it has not been explored yet
	"""
	if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		"""
		set the explored status of the current node to true
		"""
		Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		
		"""
		here we explore the out edge list of the current cluster Cx to other clusters Cy
		however, we distinguish between two types of out edges
		1) where the edge Cx->Cy was formed due to the consensus R1 relation between Cx and Cy
		these clusters Cy will go in a list of nodes "outnodes"
		2) where the edge Cx->Cy was formed due to the consensus R4 relation between Cx and Cy
		these clusters Cy were originally placed in the candidate out edge list of Cx
		here these clusters are stored in a list called "R4_outnodes"
		
		Note that individual Cy clusters should not be explored
		"""
		outnodes = []
		R4_outnodes = []
		for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
			if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				if l in Cluster_Info_Dict[root_clust_node_idx]._GetPossibleR1List():
					R4_outnodes.append(l)
				else:
					outnodes.append(l)
		
		"""
		taxa list of the current cluster
		"""
		spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		
		"""
		if R4_outnodes list is non empty, first start with an opening bracket 
		since here two or more taxa clusters will be embedded
		"""
		if (len(R4_outnodes) > 0):
			Tree_Str_List = Tree_Str_List + '('
		
		"""
		if the current cluster has no direct descendant, we just print its constituent species list
		otherwise, we enclose the newick string with additional set of opening and closing bracket
		"""
		if (len(outnodes) == 0):
			if (len(spec_list) > 1):
				Tree_Str_List = Tree_Str_List + '('
			Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			if (len(spec_list) > 1):
				Tree_Str_List = Tree_Str_List + ')'
		else:
			"""
			this is the opening bracket of additional set 
			since the cluster has one or more out clusters
			"""
			Tree_Str_List = Tree_Str_List + '('
			
			"""
			print the constituent species list
			as there is one or more out edge clusters, we do not need to check about the cardinality of the species set 
			or insert additional set of brackets
			"""
			Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			Tree_Str_List = Tree_Str_List + ','  
			
			"""
			opening bracket before printing the contents of out edge clusters of the current cluster
			"""
			Tree_Str_List = Tree_Str_List + '('
			
			"""
			navigate through individual out edge clusters connected to this cluster
			and print their contents via recursive calling of this function
			"""
			for i in range(len(outnodes)):
				if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					if (i < (len(outnodes) - 1)):
						# we check whether any subsequent node belonging to the outnodes list
						# is left for traverse
						j = i + 1
						while (j < len(outnodes)):
							if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								break
							j = j + 1	      
						# in this case, we append one comma
						if (j < len(outnodes)):
							Tree_Str_List = Tree_Str_List + ','
			
			"""
			closing bracket after printing the contents of out edge clusters of the current cluster
			"""
			Tree_Str_List = Tree_Str_List + ')'
			
			"""
			this is the closing bracket of additional set 
			since the cluster has one or more out clusters
			"""
			Tree_Str_List = Tree_Str_List + ')'
		
		"""
		if R4_outnodes list is non empty, print the contents of clusters 
		which are connected by out edges to the current cluster
		Note: all these out edge connections were originally initiated by R4 consensus relation
		"""
		if (len(R4_outnodes) > 0):
			
			"""
			opening bracket before printing the contents of R4 out edge clusters of the current cluster
			"""
			Tree_Str_List = Tree_Str_List + '('
			
			"""
			navigate through individual out edge clusters connected to this cluster
			and print their contents via recursive calling of this function
			"""
			for i in range(len(R4_outnodes)):
				if (Cluster_Info_Dict[R4_outnodes[i]]._GetExploredStatus() == 0):  
					Tree_Str_List = Tree_Str_List + PrintNewick(R4_outnodes[i])
					if (i < (len(R4_outnodes) - 1)):
						# we check whether any subsequent node belonging to the R4_outnodes list
						# is left for traverse
						j = i + 1
						while (j < len(R4_outnodes)):
							if (Cluster_Info_Dict[R4_outnodes[j]]._GetExploredStatus() == 0):  
								break
							j = j + 1
						# in this case, we append one comma
						if (j < len(R4_outnodes)):
							Tree_Str_List = Tree_Str_List + ','
		
			"""
			closing bracket after printing the contents of R4 out edge clusters of the current cluster
			"""
			Tree_Str_List = Tree_Str_List + ')'
		
		"""
		if R4_outnodes list is non empty, end with a closing bracket 
		"""
		if (len(R4_outnodes) > 0):
			Tree_Str_List = Tree_Str_List + ')'
		
	return Tree_Str_List    

##--------------------------------------------------------
#"""
#this is the original function - sourya
#this function prints the tree in Newick format
#"""
#def PrintNewick(root_clust_node_idx):
	#if 0:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList()

	#Tree_Str_List = ''
	
	#"""
	#process the node provided it has not been explored yet
	#"""
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		## set the explored status of the current node to true
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		## get the out edge list of the current node which are not explored yet 
		#outnodes = []
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#outnodes.append(l)
		## comment - sourya
		#if (len(outnodes) == 0):
		## add - sourya
		##if (len(outnodes) <= 1):
			#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		#else:
			#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
			#Tree_Str_List = Tree_Str_List + ','    
			#Tree_Str_List = Tree_Str_List + '('
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					#if (i < (len(outnodes) - 1)):
						## we check whether any subsequent node belonging to the outnodes list
						## is left for traverse
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1	      
						## in this case, we append one comma
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','
			
			#Tree_Str_List = Tree_Str_List + ')'
			#Tree_Str_List = Tree_Str_List + ')'
		
	#return Tree_Str_List    

#--------------------------------------------------------
"""
this function defines relationship between a pair of nodes in a tree
the relationship is either ancestor / descendant, or siblings, or no relationship 
@parameters: 
	wt_taxa_subset: If True, then the intersection between 
									the und_tax_list and curr_tree_taxa is accounted
	xl_val: excess gene count (normalized) between this couplet
	lca_level: level of the LCA node of this couplet
	curr_tree_taxa: set of taxa belonging to the current gene tree
"""
def DefineLeafPairReln(xl_val, lca_level, node1, node2, reln_type, curr_tree_taxa, wt_taxa_subset):

	node1_level = node1.level()
	node2_level = node2.level()

	key1 = (node1.taxon.label, node2.taxon.label)
	key2 = (node2.taxon.label, node1.taxon.label)
	
	"""
	count of taxa of the current tree
	"""
	Curr_tree_taxa_count = len(curr_tree_taxa)

	"""
	check if key2 is already in the couplet dictionary
	here, process the couplet statistics, and return
	"""
	if key2 in TaxaPair_Reln_Dict:
		if (wt_taxa_subset == True):
			"""
			here the intersect_ratio is computed by taking the intersection 
			of the nodes underlying the LCA of the couplet
			"""
			und_tax_list = TaxaPair_Reln_Dict[key2]._GetUnderlyingTaxonList()
			len_und_tax_list = len(und_tax_list)
			intersect_ratio = (len(set(und_tax_list) & set(curr_tree_taxa)) * 1.0) / len_und_tax_list
		else:
			"""
			otherwise, the ratio is computed by normalizing the current tree taxa 
			count with the total number of taxa (for all the input trees)
			"""
			intersect_ratio = (len(curr_tree_taxa) * 1.0) / len(COMPLETE_INPUT_TAXA_LIST)
		"""
		fill the entries for this couplet
		"""
		TaxaPair_Reln_Dict[key2]._AddSupportingTree()
		if (wt_taxa_subset == True):
			TaxaPair_Reln_Dict[key2]._AddXLVal(xl_val / intersect_ratio)
		else:
			TaxaPair_Reln_Dict[key2]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key2]._AddEdgeCount(Complementary_Reln(reln_type), intersect_ratio)
		#-----------------------
		if (node1_level < node2_level):
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(1, \
				((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
		elif (node1_level > node2_level):
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(0, \
				((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
		else:	#if (node1_level == node2_level):
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(2, 0)
		
		if (reln_type == RELATION_R4):
			if ((node1_level - lca_level) == 2) and ((node2_level - lca_level) > 2):
				TaxaPair_Reln_Dict[key2]._AddFreqPseudoR1(1, 1)
			elif ((node1_level - lca_level) > 2) and ((node2_level - lca_level) == 2):
				TaxaPair_Reln_Dict[key2]._AddFreqPseudoR1(0, 1)
		#-----------------------
		return

	"""
	otherwise, create the couplet dictionary key1 if it is not existing
	"""
	if key1 not in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		
	"""
	fill the couplet statistics for the entry key1
	"""
	if (wt_taxa_subset == True):
		"""
		here the intersect_ratio is computed by taking the intersection 
		of the nodes underlying the LCA of the couplet
		"""
		und_tax_list = TaxaPair_Reln_Dict[key1]._GetUnderlyingTaxonList()
		len_und_tax_list = len(und_tax_list)
		intersect_ratio = (len(set(und_tax_list) & set(curr_tree_taxa)) * 1.0) / len_und_tax_list
	else:
		"""
		otherwise, the ratio is computed by normalizing the current tree taxa 
		count with the total number of taxa (for all the input trees)
		"""
		intersect_ratio = (len(curr_tree_taxa) * 1.0) / len(COMPLETE_INPUT_TAXA_LIST)
	"""
	fill the entries for this couplet
	"""
	TaxaPair_Reln_Dict[key1]._AddSupportingTree()
	if (wt_taxa_subset == True):
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val / intersect_ratio)
	else:
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
	TaxaPair_Reln_Dict[key1]._AddEdgeCount(reln_type, intersect_ratio)
	#------------------------------------------------
	if (node1_level < node2_level):
		TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, \
			((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
	elif (node1_level > node2_level):
		TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, \
			((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
	else:	#if (node1_level == node2_level):
		TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(2, 0)

	if (reln_type == RELATION_R4):
		if ((node1_level - lca_level) == 2) and ((node2_level - lca_level) > 2):
			TaxaPair_Reln_Dict[key1]._AddFreqPseudoR1(0, 1)
		elif ((node1_level - lca_level) > 2) and ((node2_level - lca_level) == 2):
			TaxaPair_Reln_Dict[key1]._AddFreqPseudoR1(1, 1)
	#------------------------------------------------
	return

#--------------------------------------------------------
"""
this function derives couplet relations belonging to one tree
that is provided as an input argument to this function
@parameters:  
	WEIGHT_TAXA_SUBSET: If True, the relation takes care of the 
											set of taxa underlying the LCA node for this couplet
"""
def DeriveCoupletRelations(Curr_tree, WEIGHT_TAXA_SUBSET):
	"""
	taxa set of the current tree, and also the count of taxa
	"""
	curr_tree_taxa = Curr_tree.infer_taxa().labels()
	Curr_tree_taxa_count = len(curr_tree_taxa)  

	"""
	traverse the internal nodes of the tree in postorder fashion
	"""
	for curr_node in Curr_tree.postorder_internal_node_iter():
		"""
		compute the excess gene count value 
		associated with this node    
		"""
		if (WEIGHT_TAXA_SUBSET == True):
			xl_val = (len(curr_node.leaf_nodes()) - 2)
		else:
			"""
			normalized value of excess gene count 
			with respect to the number of taxa of the current tree
			"""
			xl_val = ((len(curr_node.leaf_nodes()) - 2) * 1.0) / Curr_tree_taxa_count
		
		# this is the level value associated with this node
		curr_node_level = curr_node.level()
		
		"""
		list the leaf and internal children of the current node
		"""
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		"""
		pair of leaf nodes will be related by sibling relations
		"""
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					DefineLeafPairReln(xl_val, curr_node_level, curr_node_child_leaf_nodes[i], \
						curr_node_child_leaf_nodes[j], RELATION_R3, curr_tree_taxa, WEIGHT_TAXA_SUBSET)
		
		"""
		one leaf node (direct descendant) and another leaf node (under one internal node)
		will be related by ancestor / descendant relations
		"""
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						DefineLeafPairReln(xl_val, curr_node_level, p, r, RELATION_R1, curr_tree_taxa, WEIGHT_TAXA_SUBSET)
		
		"""
		finally a pair of leaf nodes which are descendant 
		of internal nodes will be related by RELATION_R4 relation
		"""
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							DefineLeafPairReln(xl_val, curr_node_level, p, q, RELATION_R4, curr_tree_taxa, WEIGHT_TAXA_SUBSET)

#--------------------------------------------------------
"""
For a particular couplet, this function checks all the taxa 
belonging under its LCA node 
for all of the input trees supporting this couplet
"""
def FindCoupletUnderlyingTaxon(Curr_tree):
	"""
	traverse the internal nodes of the tree in postorder fashion
	check all the couplets whose LCA node is the current internal node 
	"""
	for curr_node in Curr_tree.postorder_internal_node_iter():
		"""
		taxa set belonging under this internal node
		"""
		taxa_under_curr_node = GetTaxaUnderInternalNode(curr_node)
		
		"""
		list all the leaf and internal nodes under curr_node 
		"""
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		"""
		pair of leaf nodes under curr_node
		"""
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				node1 = curr_node_child_leaf_nodes[i]
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					node2 = curr_node_child_leaf_nodes[j]
					key1 = (node1.taxon.label, node2.taxon.label)
					key2 = (node2.taxon.label, node1.taxon.label)
					if key1 in TaxaPair_Reln_Dict:
						TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
					elif key2 in TaxaPair_Reln_Dict:
						TaxaPair_Reln_Dict[key2]._AppendUnderlyingTaxonList(taxa_under_curr_node)
					else:
						"""
						establish the couplet key first
						"""
						TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
						TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
		
		"""
		one leaf node (direct descendant) and another leaf node (under one internal node)
		"""
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				node1 = p
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						node2 = r
						key1 = (node1.taxon.label, node2.taxon.label)
						key2 = (node2.taxon.label, node1.taxon.label)
						if key1 in TaxaPair_Reln_Dict:
							TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
						elif key2 in TaxaPair_Reln_Dict:
							TaxaPair_Reln_Dict[key2]._AppendUnderlyingTaxonList(taxa_under_curr_node)
						else:
							"""
							establish the couplet key first
							"""
							TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
							TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)    
		
		"""
		a pair of leaf nodes which are descendant of internal nodes 
		"""
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						node1 = p
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							node2 = q
							key1 = (node1.taxon.label, node2.taxon.label)
							key2 = (node2.taxon.label, node1.taxon.label)
							if key1 in TaxaPair_Reln_Dict:
								TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)
							elif key2 in TaxaPair_Reln_Dict:
								TaxaPair_Reln_Dict[key2]._AppendUnderlyingTaxonList(taxa_under_curr_node)
							else:
								"""
								establish the couplet key first
								"""
								TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
								TaxaPair_Reln_Dict[key1]._AppendUnderlyingTaxonList(taxa_under_curr_node)

	return

#--------------------------------------------------------
#""" 
#this function prints the elements of the queue (which stores the couplet scores 
#for individual relations 
#"""
#def PrintQueueInfo(inp_queue, Output_Text_File):
	#fp = open(Output_Text_File, 'a')
	#for elem in inp_queue:
		#fp.write(' ' + str(elem))
	#fp.close()

#-----------------------------------------------------
"""
this function reads the input tree list file
@parameters: 
	ROOTED_TREE - whether the treelist to be read as rooted format
	PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
	INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
	INPUT_FILENAME: file containing the input treelist
"""
def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, \
		schema=INPUT_FILE_FORMAT, preserve_underscores=PRESERVE_UNDERSCORE, \
			default_as_rooted=ROOTED_TREE)

	return Inp_TreeList

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

#-----------------------------------------------------
# this is the taxa list generated from current internal node
def GetTaxaUnderInternalNode(curr_node):
	taxa_list_from_curr_internal_node = []
	for n in curr_node.leaf_nodes():
		taxa_list_from_curr_internal_node.append(n.taxon.label)
	return taxa_list_from_curr_internal_node

#----------------------------------------
def Complementary_Reln(inp_reln):
	if (inp_reln == RELATION_R3) or (inp_reln == RELATION_R4):
		return inp_reln
	elif (inp_reln == RELATION_R1):
		return RELATION_R2
	else:
		return RELATION_R1

#------------------------------------------------
"""
this function computes average XL information between a pair of taxa clusters
@param: taxa_clust1: first taxa list
				taxa_clust2: second taxa list
				DIST_MAT_TYPE: Type of distance employed
				single_elem: can contain one of possible three values
				0: only one element of taxa_clust1 and one element of taxa_clust2 will be compared
				1: cluster containing taxa_clust1[0] and cluster containing taxa_clust2[0] will be compared
				2: All pairs of elements of taxa_clust1 and taxa_clust2 will be compared
"""
def FindAvgXL(taxa_clust1, taxa_clust2, DIST_MAT_TYPE, single_elem=2, type_of_output=0):
	"""
	if single_elem = 0
	we compare taxa_clust1[0] and taxa_clust2[0], in terms of the preorder level
	
	if single_elem = 1
	we check the first preorder level taxon of both lists taxa_clust1 and taxa_clust2
	suppose the taxon names are taxa1 and taxa2
	but instead of comparing taxa1 and taxa2 only
	we compare the original taxa clusters (may have cardinality > 1) containing taxa1 and taxa2
	
	if single_elem = 2
	we compare pairwise all the elements belonging to taxa_clust1 and taxa_clust2
	"""
	if (single_elem == 1):
		taxa1 = taxa_clust1[0]
		taxa2 = taxa_clust2[0]
		clust1 = Taxa_Info_Dict[taxa1]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[taxa2]._Get_Taxa_Part_Clust_Idx()
		taxa_list1 = Cluster_Info_Dict[clust1]._GetSpeciesList()
		taxa_list2 = Cluster_Info_Dict[clust2]._GetSpeciesList()
	elif (single_elem == 2):
		taxa_list1 = taxa_clust1
		taxa_list2 = taxa_clust2 
	else:
		taxa_list1 = []
		taxa_list1.append(taxa_clust1[0])
		taxa_list2 = []
		taxa_list2.append(taxa_clust2[0])
		
	curr_taxa_pair_list = []
	for x1 in taxa_list1:
		for x2 in taxa_list2:  
			key1 = (x1, x2)
			key2 = (x2, x1)
			#print 'key1: ', key1, ' key2: ', key2
			if key1 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				curr_taxa_pair_list.append(val)
			elif key2 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				curr_taxa_pair_list.append(val)
	
	# average of this pairwise list is used as the XL approximation
	if (len(curr_taxa_pair_list) > 0):
		if (type_of_output == 0):
			return (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
		else:
			return max(curr_taxa_pair_list)
			#return min(curr_taxa_pair_list)
	else:
		return 0
	
#-----------------------------------------------------------------
"""
this function returns the frequency of R1 relation from taxa1 to taxa2 (or frequency of R2 relation from taxa2 to taxa1)
@param: if percent_tree_frac is 1, the output is normalized with the number of supporting trees
"""
def GetR1Freq(taxa1, taxa2, percent_tree_frac=False):
	val = 0
	key1 = (taxa1, taxa2)
	key2 = (taxa2, taxa1)
	if key1 in TaxaPair_Reln_Dict:
		val = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
		if (percent_tree_frac == True):
			val = (val * 1.0) / TaxaPair_Reln_Dict[key1]._GetNoSupportTrees()
	elif key2 in TaxaPair_Reln_Dict:
		val = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
		if (percent_tree_frac == True):
			val = (val * 1.0) / TaxaPair_Reln_Dict[key2]._GetNoSupportTrees()
	
	return val

#-----------------------------------------------------------------
"""
this function returns the list of taxa underlying the given internal node
in preorder traversal
@param: inp_node: Input node under which the taxa set will be explored
				taxa_list: Output taxa list in preorder traversal order
				inp_set_of_taxa: A superset of taxon; the 'taxa_list' should be a subset of it
"""
def GetPreorderTaxaList(inp_node, taxa_list, inp_set_of_taxa):
	for n in inp_node.preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in inp_set_of_taxa:
				taxa_list.append(n.taxon.label)
	
	return taxa_list
