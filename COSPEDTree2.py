#!/usr/bin/env python


##---------------------------------------------
""" 
this program is used to generate a supertree (consensus) from a set of constituent trees
the input is multiple source trees
each of the trees need to be decomposed to couplets 
then the couplets need to be joined
there may be conflicts among the input tree - we have to select the consensus

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 15.01.2014 - public release
V2.0 - 15.03.2014 - rewritten for lowering time complexity
V3.0 - 28.02.2015 - added binary suprtree and level based scoring options
V4.0 - 30.03.2015 - added fast input tree processing
V5.0 - 17.06.2015 - binary refinement and github release
V6.0 - 01.03.2016 - modified the code to incorporate cluster based connections
""" 

# Copyright 2013, 2014, 2015, 2016 
# Sourya Bhattacharyya and Jayanta Mukherjee.
# All rights reserved.
#
# See "LICENSE.txt" for terms and conditions of usage.
#
#---------------------------------------------

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import Process_Queues
from Process_Queues import *
import ReachGraph_Update
from ReachGraph_Update import *
import UtilFunc
from UtilFunc import *
import RefineTree
from RefineTree import *
import Conflict_Detect
from Conflict_Detect import *
import Cluster_Manage
from Cluster_Manage import *

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
	parser = OptionParser()

	parser.add_option("-I", "--INPFILE", \
				type="string", \
				action="store", \
				dest="INP_FILENAME", \
				default="", \
				help="name of the input file (along with complete or relative path) containing input trees")

	parser.add_option("-O", "--OUTFILE", \
				type="string", \
				action="store", \
				dest="OUT_FILENAME", \
				default="", \
				help="name of the output file which will contain the output supertree")  

	parser.add_option("-p", "--inpform", \
				type="int", \
				action="store", \
				dest="inp_file_format", \
				default=1, \
				help="1 - input file format is NEWICK (default) \
				2 - input file format is NEXUS")

	parser.add_option("-b", "--binary", \
				action="store_true", \
				dest="binary_suptr", \
				default=False, \
				help="if TRUE, it produces a strictly binary supertree. \
				Otherwise, the tree can be non-binary. Default FALSE.")
					
	parser.add_option("-u", "--underscore", \
				action="store_false", \
				dest="preserve_underscores", \
				default=True, \
				help="this is a boolean flag option \
				using this option toggles the existing configuration (Default TRUE) \
				if TRUE, then this option preserves the underscores of the names of taxa \
				so, enabling this option do not preserve the underscores")  
					
	#parser.add_option("-n", "--njrule", \
				#type="int", \
				#action="store", \
				#dest="NJ_type", \
				#default=2, \
				#help="valid only if binary supertree is produced \
				#1 - classical NJ method \
				#2 - Normalized couplet statistic for agglomeration (Default)")     

	parser.add_option("-n", "--njmerge", \
				type="int", \
				action="store", \
				dest="nj_merge_clust", \
				default=1, \
				help="1 - Use earlier cluster merge\
				2 - Use new cluster merge")     

	parser.add_option("-w", "--weighttaxa", \
				action="store_false", \
				dest="weight_taxa_subset", \
				default=True, \
				help="this is a boolean flag option \
				using this option toggles the existing configuration (Default TRUE) \
				if TRUE, then this option weighs couplet statistics according \
				to the size of taxa subset underlying MRCA of that couplet")  

	#parser.add_option("-x", "--mppsolvemetric", \
				#type="int", \
				#action="store", \
				#dest="MPP_solve_metric", \
				#default=2, \
				#help="1 - Use priority measure based selection (higher value)\
				#2 - Use XL based selection (lower value)")     

	parser.add_option("-d", "--distmat", \
				type="int", \
				action="store", \
				dest="dist_mat_type", \
				default=1, \
				help="1 - Mean of XL \
				2 - Mean(Average, Mode based Avg) of XL")     

	opts, args = parser.parse_args()
	return opts, args
  
  
#-----------------------------------------------------
"""
main function of the source code
"""
def main():  
	opts, args = parse_options()

	ROOTED_TREE = False 	#True #opts.default_rooted
	PRESERVE_UNDERSCORE = opts.preserve_underscores
	if (opts.inp_file_format == 1):
		INPUT_FILE_FORMAT = 'newick'
	else:
		INPUT_FILE_FORMAT = 'nexus'
	INPUT_FILENAME = opts.INP_FILENAME
	OUTPUT_FILENAME = opts.OUT_FILENAME
	BINARY_SUPERTREE_OPTION = opts.binary_suptr
	NJ_RULE_USED = 2	#opts.NJ_type 
	WEIGHT_TAXA_SUBSET = opts.weight_taxa_subset
	MPP_SOLVE_METRIC = 2	#opts.MPP_solve_metric
	DIST_MAT_TYPE = opts.dist_mat_type
	NJ_MERGE_CLUST = 1	#2	#opts.nj_merge_clust
	NO_OF_QUEUES = 1	# 2

	if (INPUT_FILENAME == ""):
		print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
		return
	else:
		print 'input filename: ', INPUT_FILENAME
	
	"""
	according to the location of input filename
	adjust the locations of the output files as well
	"""
	k = INPUT_FILENAME.rfind("/")
	if (k == -1):
		dir_of_inp_file = './'
	else:
		dir_of_inp_file = INPUT_FILENAME[:(k+1)]
		
	input_file_name = INPUT_FILENAME[(k+1):]
	if (DEBUG_LEVEL > 1):
		print 'dir_of_inp_file: ', dir_of_inp_file  
		
	#-------------------------
	if (OUTPUT_FILENAME == ""):
		"""
		first create the output directory containing the results
		the directory name should also reflect the settings of input parameters
		"""
		dir_of_curr_exec = dir_of_inp_file + 'COSPEDBTree'
		
		"""
		according to the settings of input parameters, 
		output directory name is customized
		"""
		if (BINARY_SUPERTREE_OPTION == True):
			dir_of_curr_exec = dir_of_curr_exec + '_B_1'
		else:
			dir_of_curr_exec = dir_of_curr_exec + '_B_0'
		if (WEIGHT_TAXA_SUBSET == True):
			dir_of_curr_exec = dir_of_curr_exec + '_W_1'
		else:
			dir_of_curr_exec = dir_of_curr_exec + '_W_0'
		if (BINARY_SUPERTREE_OPTION == True):
			dir_of_curr_exec = dir_of_curr_exec + '_N_' + str(NJ_RULE_USED) + '_D_' + str(DIST_MAT_TYPE)
		""" 
		create the directory
		"""
		if (os.path.isdir(dir_of_curr_exec) == False):
			mkdr_cmd = 'mkdir ' + dir_of_curr_exec
			os.system(mkdr_cmd)               
		
	else:
		"""
		here the output filename is specified
		so extract the output directory from that file name
		"""
		k = OUTPUT_FILENAME.rfind("/")
		if (k == -1):
			dir_of_curr_exec = './'
		else:
			dir_of_curr_exec = OUTPUT_FILENAME[:(k+1)]

	"""
	decide the output text file containing the text results
	"""
	Output_Text_File = dir_of_curr_exec + '/' + 'COSPEDBTree_Complete_Desription.txt'
	if (DEBUG_LEVEL > 1):
		print 'dir_of_curr_exec: ', dir_of_curr_exec  
		print 'Output_Text_File: ', Output_Text_File      

	#-------------------------
	# open the output text file
	fp = open(Output_Text_File, 'w')    
	
	fp.write('\n ================ status of options ================= (1 means ON)')
	fp.write('\n ROOTED_TREE: ' + str(ROOTED_TREE))
	fp.write('\n PRESERVE_UNDERSCORE: ' + str(PRESERVE_UNDERSCORE))
	fp.write('\n BINARY SUPERTREE OPTION: ' + str(BINARY_SUPERTREE_OPTION))
	fp.write('\n ===>>>  processing the file now ======== ')

	# note the program beginning time 
	start_timestamp = time.time()
		
	#-------------------------------------  
	""" 
	read the source trees collection and store it in a treelist
	individual elements of this collection is a source tree 
	"""
	Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  

	#-------------------------------------  
	""" 
	from the input trees, note the number of taxa (total)
	and also define the class instances corresponding to single taxa
	"""
	for tr_idx in range(len(Input_Treelist)):
		taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
		if (DEBUG_LEVEL > 1):
			fp.write('\n Tree no : ' + str(tr_idx+1) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
		if (DEBUG_LEVEL > 2):
			fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
		for i in range(len(taxa_labels_curr_tree)):
			if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
				COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])

	"""
	we also define one structure "Taxa_Info_Dict" marked by a taxa
	this defines a single taxon, and becomes an instance of the class "Single_Taxa"
	"""
	for label in COMPLETE_INPUT_TAXA_LIST:
		Taxa_Info_Dict.setdefault(label, Single_Taxa())

	"""
	total number of taxa within the input tree set
	"""
	number_of_taxa = len(COMPLETE_INPUT_TAXA_LIST)  

	if (DEBUG_LEVEL >= 0):
		fp.write('\n  total no of taxa: ' + str(number_of_taxa))
		fp.write('\n  total no of trees: ' + str(len(Input_Treelist)))

	data_read_timestamp1 = time.time()	# note the timestamp

	if (DEBUG_LEVEL >= 0):
		fp.write('\n Time taken to read the taxa information: ' \
			+ str(data_read_timestamp1 - start_timestamp))
		
	#---------------------------
	"""
	for individual couplets, this function lists all the taxa 
	lying within their LCA node for all the input trees
	"""
	if (WEIGHT_TAXA_SUBSET == True):
		for tr_idx in range(len(Input_Treelist)):
			FindCoupletUnderlyingTaxon(Input_Treelist[tr_idx])
			
	data_read_timestamp2 = time.time()	# note the timestamp
	
	if (DEBUG_LEVEL >= 0):
		fp.write('\n Time taken to compute the underlying taxon for weight \
			computation: ' + str(data_read_timestamp2 - data_read_timestamp1))
	#---------------------------
	"""
	now process individual trees to find the couplet relations within those trees
	"""
	for tr_idx in range(len(Input_Treelist)):
		DeriveCoupletRelations(Input_Treelist[tr_idx], WEIGHT_TAXA_SUBSET)

	if (DEBUG_LEVEL > 2):
		fp.write('\n len Taxa_Info_Dict: ' + str(len(Taxa_Info_Dict)))
		fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
		fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
	
	# note the timestamp
	data_read_timestamp = time.time()
		
	if (DEBUG_LEVEL >= 0):
		fp.write('\n Time taken to read the couplet relations: ' + str(data_read_timestamp - data_read_timestamp2))

	# close the output text file
	fp.close()
	
	#------------------------------------------------------------
	""" 
	Here we allocate the list / storage space of taxa clusters
	each taxa cluster is supposed to store the taxa subsets 
	related via relation r3 (simultaneous speciation)
	initially all the clusters contain one taxa
	each of the cluster has the index of the corresponding 
	taxa in the COMPLETE_INPUT_TAXA_LIST 
	"""
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		Create_Cluster_Taxa_Label(i, COMPLETE_INPUT_TAXA_LIST[i])
	
	#------------------------------------------------------------
	""" 
	we initialize the Reachability_Graph_Mat
	dimension: N X N where N = number of taxa clusters = number of taxa (initially)
	this is a numpy 2D array 
	values Mat[x][y] = 1 means x->y
	Mat[x][y] = Mat[y][x] = 2 means x and y are connected via relation r4
	as taxa subsets are clustered, no of distinct taxa clusters reduces
	so, this matrix dimension is also decreased
	"""
	Reachability_Graph_Mat = numpy.zeros((number_of_taxa, number_of_taxa), dtype=numpy.int)

	fp = open(Output_Text_File, 'a')
	if (DEBUG_LEVEL > 0):
		fp.write('\n shape of Reachability_Graph_Mat: ' \
			+ str(numpy.shape(Reachability_Graph_Mat)))

	# this information is printed to know the maximum possible 
	# iterations that the while loops will undergo
	if (DEBUG_LEVEL > 1):
		fp.write('\n length of conflicting priority queue : ' \
			+ str(len(Cost_List_Taxa_Pair_Multi_Reln)))

	fp.close()
	#------------------------------------------------------------
	""" 
	here, we process all the couplets to extract couplet statistics
	"""  
	for l in TaxaPair_Reln_Dict:
		"""
		single_edge_exist: if TRUE, means that only one type of relation 
		is supported (with respect to input trees) between this couplet
		detection of it during setting the priority values of different relations
		basically, we are looking for the consensus relation
		"""
		single_edge_exist_list = TaxaPair_Reln_Dict[l]._CheckNonConflictingCouplet()
		single_edge_exist = single_edge_exist_list[0]
		consensus_reln_type = single_edge_exist_list[1]

		"""
		first set the priority values of this couplet
		"""
		TaxaPair_Reln_Dict[l]._SetConnPrVal()

		""" 
		calculate the support score and priority measures for individual couplets
		and for individual relations
		it is the product of frequency and the priority measures 
		"""
		TaxaPair_Reln_Dict[l]._SetCostMetric()
		
		#"""
		#sort the frequency list of this couplet, in descending order
		#"""
		#TaxaPair_Reln_Dict[l]._SortFreqList()
		
		if (NO_OF_QUEUES == 2) and (single_edge_exist == True):
			"""
			if we are allowed to use two different queues for storing conflicting and non-conflicting couplets
			
			for this particular couplet, we select only the consensus relation 
			and a relation having frequency very close to the consensus relation
			to insert in the support score queue
			"""
			"""
			add this relation as the allowed relation for this couplet
			"""
			TaxaPair_Reln_Dict[l]._AddAllowedReln(consensus_reln_type)
			"""
			support score measure for this particular relation
			if support_score > 0, then the relation is majority consensus
			"""
			support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(consensus_reln_type)
			"""
			add this score in the queue
			Note: we have only one queue for processing support scores
			"""
			sublist = [l[0], l[1], consensus_reln_type, support_score]
			Cost_List_Taxa_Pair_Single_Reln.append(sublist)
		else:
			"""
			for a particular couplet, we select only the consensus relation 
			and a relation having frequency very close to the consensus relation
			to insert in the support score queue
			"""
			#max_freq = TaxaPair_Reln_Dict[l]._GetConsensusFreq()
			for reln_type in range(4):
				curr_reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type)
				#if (TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type) >= (PERCENT_MAX_FREQ * max_freq)):
				#if (TaxaPair_Reln_Dict[l]._CheckMaxor2ndMaxReln(curr_reln_freq) == 1):
				if (curr_reln_freq > 0):
					"""
					add this relation as the allowed relation for this couplet
					"""
					TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
					"""
					support score measure for this particular relation
					if support_score > 0, then the relation is majority consensus
					"""
					support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(reln_type)
					"""
					add this score in the queue
					Note: we have only one queue for processing support scores
					"""
					sublist = [l[0], l[1], reln_type, support_score]
					Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
					
					#----------------------------------------
					#"""
					#if the relation type is R3
					#and the support score is also positive (majority consensus relation)
					#then add the relation R3 between this couplet 
					#within the list of allowed relations in the final supertree

					#Note: the relation should not be established earlier
					#or the relation should not generate any conflict
					
					#update the Reachability_Graph_Mat also
					#"""
					#if (reln_type == RELATION_R3) and (support_score > 0) and (single_edge_exist == True):
						#clust1 = Taxa_Info_Dict[l[0]]._Get_Taxa_Part_Clust_Idx()
						#clust2 = Taxa_Info_Dict[l[1]]._Get_Taxa_Part_Clust_Idx()
						#if (CheckExistingConn(clust1, clust2, Reachability_Graph_Mat, RELATION_R3, Output_Text_File) == 0):
							#if (CheckTransitiveConflict(clust1, clust2, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
								#if (DEBUG_LEVEL >= 2):
									#fp = open(Output_Text_File, 'a')    
									#fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' \
										#+ str(l[0]) + ' and ' + str(l[1]) + ' relation type: ' + str(RELATION_R3))
									#fp.close()
								## transitive update of the Reachability_Graph_Mat
								#Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, clust1, clust2, RELATION_R3)
		
	#------------------------------------------------------------
	# we print the original connection status for all the tree nodes
	if (DEBUG_LEVEL >= 2):
		for l in TaxaPair_Reln_Dict:
			#print 'printing info for the TaxaPair_Reln_Dict key: ', l
			TaxaPair_Reln_Dict[l]._PrintRelnInfo(l, Output_Text_File)

	#------------------------------------------------------------
	## add - sourya
	#"""
	#here we check whether any three consensus relations are conflicting
	#in such a case, we remove one of the consensus relations from the support score queue
	#"""
	#Check_Consensus_Reln_Conflict()
	## end add - sourya
	#------------------------------------------------------------
	""" 
	sort the priority queues according to the support score value of individual relations
	4th field of individual sublist denotes the support score
	we use max priority queue implementation for storing and sorting these queues
	"""
	if (NO_OF_QUEUES == 2) and (len(Cost_List_Taxa_Pair_Single_Reln) > 0):
		# sorting the non-conflicting queue
		Sort_Priority_Queue(Cost_List_Taxa_Pair_Single_Reln)
	
	# sorting the conflicting queue
	Sort_Priority_Queue(Cost_List_Taxa_Pair_Multi_Reln)

	data_initialize_timestamp = time.time()	# note the timestamp
	#------------------------------------------------------------
	"""
	first we process the non-conflicting queue
	provided we use it
	"""
	if (NO_OF_QUEUES == 2) and (len(Cost_List_Taxa_Pair_Single_Reln) > 0):
		Reachability_Graph_Mat = Proc_Queue_Pos_Score(Reachability_Graph_Mat, Output_Text_File, 1)
	
	"""
	we process the conflicting queue now
	"""
	Reachability_Graph_Mat = Proc_Queue_Pos_Score(Reachability_Graph_Mat, Output_Text_File, 0)
	
	#------------------------------------------------------------
	# we print the final connection status for all the tree nodes
	if (DEBUG_LEVEL > 2):
		for l in Taxa_Info_Dict:
			#print 'printing information for the Taxa ', l
			Taxa_Info_Dict[l]._PrintFinalTaxaInfo(l, Output_Text_File) 

	# print the cluster information 
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
		fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
		fp.write(str(CURRENT_CLUST_IDX_LIST))
		fp.write('\n ========== cluster information after reachability graph generation =============')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)

	# note the timestamp
	reachability_graph_form_timestamp = time.time()  
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n Contents of Candidate_Out_Edge_Cluster_List: ' \
			+ str(Candidate_Out_Edge_Cluster_List)) 
		fp.close()
	#------------------------------------------------------------
	"""
	now we process all the elements of the Candidate_Out_Edge_Cluster_List one by one
	"""
	Reachability_Graph_Mat = Process_Candidate_Out_Edge_Cluster_List(Reachability_Graph_Mat, DIST_MAT_TYPE, Output_Text_File)
	
	# print the cluster information 
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
		fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
		fp.write(str(CURRENT_CLUST_IDX_LIST))
		fp.write('\n ========== cluster information after processing Candidate_Out_Edge_Cluster_List =============')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
	#------------------------------------------------------------
	""" 
	now perform the transitive reduction of the closure formed by 
	connection of the cluster of nodes in the above operation
	this is required to handle the following scenario:
	suppose, there exists a case such that A->C, B->C and A->B
	then in the final graph, only A->B and B->C information needs to be preserved
	in order to form the DAG 
	"""
	CompressDirectedGraph(Reachability_Graph_Mat)
	
	# print the cluster information 
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
		fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
		fp.write(str(CURRENT_CLUST_IDX_LIST))
		fp.write('\n ========== cluster information after transitive reduction =============')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
	
	#------------------------------------------------------------
	# now instead of arbitrary assignment of the parent node for individual clusters 
	# we assign parent node according to the source tree relationships
	# this will solve the multiple parent problem C2 as discussed in the manuscript 
	# this is a new addition and marked under the DFS based parent refinement option 
	SolveMultipleParentC2Problem(Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE)
				
	# print the cluster information 
	if (DEBUG_LEVEL > 0):
		fp = open(Output_Text_File, 'a')
		fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
		fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
		fp.write(str(CURRENT_CLUST_IDX_LIST))    
		fp.write('\n ========== cluster information after solving multiple parent problem =============')
		fp.close()
		for i in Cluster_Info_Dict:
			#print 'printing the information for cluster node: ', i
			Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
			
	# note the timestamp
	cluster_of_node_refine_species_timestamp1 = time.time()  

	#----------------------------------------------
	""" 
	construct supertree from the generated DAG 
	scheme: repeatedly extract the nodes (taxa clusters) with minimum indegree
	print taxa information constituent within this cluster
	maintain appropriate paranthesis according to the newick format
	"""
	""" 
	the variable no_of_components is used in association with the used depth first technique 
	if no_of_components > 1, it signifies that a forest has been created
	in other words, the problem C3 (no parent problem exists)
	"""
	no_of_components = 0
	while (1):
		root_clust_node_idx = Extract_Node_Min_Indeg(len(CURRENT_CLUST_IDX_LIST))
		if (root_clust_node_idx == -1):
			break
		Tree_Str = PrintNewick(root_clust_node_idx)
		no_of_components = no_of_components + 1
		if (no_of_components == 1):	# first component
			Final_Supertree_Str = Tree_Str
		else:
			Final_Supertree_Str = Final_Supertree_Str + ',' + Tree_Str

	""" 
	with the final tree string, finally generate the tree result 
	this is also required to tackle the creation of a forest during the DFS (no parent problem)
	this procedure adheres to the basic COSPEDTree mechanism
	"""
	Final_Supertree_Str = '(' + Final_Supertree_Str + ')'
	#---------------------------------------------------------
	fp = open(Output_Text_File, 'a')
	fp.write('\n\n\n\n **** original supertree as newick string --- ' + Final_Supertree_Str) 

	## comment - sourya
	#Final_Supertree_Str = Remove_Extra_Paranthesis(Final_Supertree_Str)  
	#fp.write('\n --- after removing extra paranthesis -- supertree as newick string --- ' + Final_Supertree_Str) 
	#fp.close()
	## end comment - sourya

	# now read this super string in a supertree containing all the input taxa
	Supertree_Final = dendropy.Tree.get_from_string(Final_Supertree_Str, \
		schema="newick")	#preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=ROOTED_TREE)
	if 0:
		Supertree_Final.print_plot()  
		
	# note the timestamp
	newick_str_formation_timestamp = time.time()  

	fp = open(Output_Text_File, 'a')
	fp.write('\n initial supertree: ' + Supertree_Final.as_newick_string())    
	Supertree_Final.update_splits(delete_outdegree_one=True)
	fp.write('\n\n\n ***** after update splits: supertree: ' + Supertree_Final.as_newick_string())    
	fp.close()
	#-------------------------------------------------------------
	if (BINARY_SUPERTREE_OPTION == True):    
		""" 
		this function removes all multifurcating clusters and produces binary tree 
		it also solves the problem C3, as mentioned in the manuscript
		"""
		Refine_Supertree_Binary_Form(Supertree_Final, Output_Text_File, \
			NJ_RULE_USED, DIST_MAT_TYPE, NJ_MERGE_CLUST)
		
		fp = open(Output_Text_File, 'a')
		fp.write('\n --- user provided option for producing strict binary supertree')
		fp.write('\n --- after binary refinement --- output tree without branch length (in newick format): ' \
			+ Supertree_Final.as_newick_string())    
		fp.close()
	else:
		fp = open(Output_Text_File, 'a')
		fp.write('\n --- user did not provide option for producing strict binary supertree - so output tree can be non-binary')
		fp.close()
		
	# write this tree on a separate text file
	if (OUTPUT_FILENAME == ""):
		out_treefilename = dir_of_curr_exec + '/' + 'cospedbtree_newick.tre'
	else:
		out_treefilename = OUTPUT_FILENAME

	Supertree_Final.write_to_path(out_treefilename, 'newick')
	
	#---------------------------------
	"""
	Note: This seems to be a redundant operation, but an important one
	Somehow, the structure "Supertree_Final" fails to read the underscores properly
	so we read the tree from the file itself, and store it again in the "Supertree_Final" structure
	"""
	Supertree_Final = dendropy.Tree.get_from_path(out_treefilename, schema='newick', \
		preserve_underscores=PRESERVE_UNDERSCORE)

	# final timestamp
	data_process_timestamp = time.time()      
	#---------------------------------
	timing_text_file = dir_of_curr_exec + '/' + 'timing_info.txt'
	fp = open(timing_text_file, 'w')  
	fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds) ')
	fp.write('\n \n reading the data: ' + str(data_read_timestamp - start_timestamp) + \
	'\n initialization of the structure: ' + str(data_initialize_timestamp - data_read_timestamp) + \
	'\n formation of the reachability graph (cluster) (after loop): ' + \
				str(reachability_graph_form_timestamp - data_initialize_timestamp) + \
	'\n multiple parent (related) problem: ' + \
				str(cluster_of_node_refine_species_timestamp1 - reachability_graph_form_timestamp) + \
	'\n newick string formation: ' + str(newick_str_formation_timestamp - cluster_of_node_refine_species_timestamp1) + \
		'\n binary tree construction: ' + str(data_process_timestamp - newick_str_formation_timestamp))
	fp.write('\n \n Total time taken (in seconds) : ' + str(data_process_timestamp - start_timestamp))  
	fp.close()

	#----------------------------------------------
	# Performance metric code
	#----------------------------------------------
	# open the output text file
	outtextfile = dir_of_curr_exec + '/' + 'FP_FN_RF_Perf.txt'
	fp = open(outtextfile, 'w')

	"""
	examine each of the source trees and find the FP, FN and RF distance 
	with respect to the generated supertree  
	"""
	sumFP = sumFN = sumRF = 0  
	sumLenSrcTree = 0
	sum_symmetric_diff = 0
	fp.write('\n \n\n total edges of supertree: ' + str(len(Supertree_Final.get_edge_set())))  

	#print 'taxon set of supertree: ', Supertree_Final.infer_taxa()
	for inp_tree_idx in range(len(Input_Treelist)):
		Curr_src_tree = Input_Treelist[inp_tree_idx]
		curr_src_tree_taxa = Curr_src_tree.infer_taxa().labels()
		curr_src_tree_no_of_taxa = len(curr_src_tree_taxa)
		
		# according to the taxa set of the current source tree, 
		# prune the supertree to get the tree portion containing only this taxa set
		pruned_tree = dendropy.Tree(Supertree_Final)
		pruned_tree.retain_taxa_with_labels(curr_src_tree_taxa)
		
		# source tree number of edges calculation
		# it is used to compute normalized RF metric values
		lenSrcTree = len(Curr_src_tree.get_edge_set())
		sumLenSrcTree = sumLenSrcTree + lenSrcTree
		fp.write('\n src tree : ' + str(Curr_src_tree))
		fp.write('\n pruned supertree: ' + str(pruned_tree))
		fp.write('\n src tree len: ' + str(lenSrcTree) + ' pruned supertree len: ' \
			+ str(len(pruned_tree.get_edge_set())))
		
		# determine the false positives and the false negatives 
		tup = Curr_src_tree.false_positives_and_negatives(pruned_tree)
		fp.write('   FP_int: ' + str(tup[0]) + '  FN_int:  ' + str(tup[1]))
		sumFP = sumFP + tup[0]
		sumFN = sumFN + tup[1]
		sumRF = sumRF + ((tup[0] + tup[1]) / 2.0)
		
		symm_diff = Curr_src_tree.symmetric_difference(pruned_tree)
		fp.write('   Symmetric difference: ' + str(symm_diff))
		sum_symmetric_diff = sum_symmetric_diff + symm_diff

	# final normalized sumFP's are computed by dividing with the number of trees
	normsumFP = (sumFP * 1.0) / sumLenSrcTree
	normsumFN = (sumFN * 1.0) / sumLenSrcTree
	normsumRF = (sumRF * 1.0) / sumLenSrcTree
	norm_symm_diff = sum_symmetric_diff / (2.0 * sumLenSrcTree)
		
	# print the final result
	fp.write('\n\n\n ===============>>>>>>>>>>>>>>> FINAL RESULTS \n \n')
	fp.write('\n ******* absolute sumFP: ' + str(sumFP) + \
		'\n ******* absolute sumFN: ' + str(sumFN) + \
		'\n ******* absolute sumRF: ' + str(sumRF) + \
		'\n ******* absolute Symmetric difference: ' + str(sum_symmetric_diff))

	fp.write('\n ===============>>>>>>>>>>>>>>> IN TERMS OF NORMALIZED \
		(DIVIDED BY THE SUM OF INTERNAL EDGES OF THE SOURCE TREES) ''')  
	fp.write('\n normsumFP: ' + str(normsumFP) + '\n normsumFN: ' + str(normsumFN) + \
		'\n normsumRF: ' + str(normsumRF) + '\n norm Symmetric Diff: ' + str(norm_symm_diff)) 

	fp.close()
	#----------------------------------------------
	# end Performance metric code
	#----------------------------------------------
	#--------------------------------------------------------------  
	# delete the storage variables associated with the current execution 

	# clear the dictionaries
	Cluster_Info_Dict.clear()
	Taxa_Info_Dict.clear()
	TaxaPair_Reln_Dict.clear()

	# clear the lists associated
	if (len(Cost_List_Taxa_Pair_Multi_Reln) > 0):
		Cost_List_Taxa_Pair_Multi_Reln[:] = []
	if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
		COMPLETE_INPUT_TAXA_LIST[:] = []
	if (len(CURRENT_CLUST_IDX_LIST) > 0):
		CURRENT_CLUST_IDX_LIST[:] = []

	# free the reachability graph (numpy array)
	del Reachability_Graph_Mat
  
##-----------------------------------------------------

if __name__ == "__main__":
    main() 
  
