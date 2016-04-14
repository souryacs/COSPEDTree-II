#!/usr/bin/env python

import dendropy
from dendropy import Tree, Taxon, TaxonSet, Node
import numpy
import time
import os
from optparse import OptionParser
#import sys
import math

# we define custom edge types
RELATION_R3 = 0	# equality relationship
RELATION_R1 = 1
RELATION_R2 = 2
RELATION_R4 = 3	# no relationship
UNDEFINED_RELATION = 4

# variables used to denote whether we use traditional NJ method
# or use a variant of it, namely the agglomerative clustering
TRADITIONAL_NJ = 1
AGGLO_CLUST = 2

## variables used to denote the distance metric employed for agglomerative clustering
#EXTRA_TAXA = 1
#PROD_EXTRA_TAXA_BRANCH_COUNT = 2
#PROD_EXTRA_TAXA_ACC_RANK = 3

#---------------------------------------------
""" 
this is a dictionary for storing information about individual taxa clusters
each cluster is basically a collection of taxa related via relation r3
"""
Cluster_Info_Dict = dict()

""" 
the dictionary defines one particular taxa and its associated information
"""
Taxa_Info_Dict = dict()

""" 
this dictionary defines the taxa pair (couplet) relations and associated operations
each entry of this dictionary is indexed by a pair of taxon labels 
"""
TaxaPair_Reln_Dict = dict()

"""
queue storing relations of non-conflicting couplets (supporting only one type of relation 
in the input trees)
"""
Cost_List_Taxa_Pair_Single_Reln = []

"""
queue storing relations of conflicting couplets (supporting more than one type of relation 
in the input trees)
"""
Cost_List_Taxa_Pair_Multi_Reln = []

""" 
this list contains the complete set of taxa present in the input trees 
"""
COMPLETE_INPUT_TAXA_LIST = []

""" 
this list contains the current set of active taxa cluster (indices)
"""
CURRENT_CLUST_IDX_LIST = []

"""
this is the debug level
set for printing the necessary information
"""
DEBUG_LEVEL = 2

# add - sourya
MODE_PERCENT = 0.25	#0.4 #0.35	#0.5 
MODE_BIN_COUNT = 40

"""
this is a threshold corresponding to the selection of R1 or R2 relation
for a non conflicting couplet with negative support score of the corresponding relation
"""
R1R2Reln_MAJ_THRS_high = 0.75	#0.8
R1R2Reln_MAJ_THRS_low = 0.65	#0.7	#0.8
R1R2Reln_MAJ_THRS_very_low = 0.6	#0.65

"""
for a couplet, relations with frequency of this percentage of the 
max (consensus) frequency, will be included in the score queue
currently we employ 50% frequency threshold
"""
PERCENT_MAX_FREQ = 0.7

"""
this list contains the set of clusters 
which need to be processed for setting at least one directed out edge to another cluster
"""
Candidate_Out_Edge_Cluster_List = []

##-----------------------------------------------------
""" 
this class defines a leaf node of input candidate source trees
that is, it corresponds to one taxa 
"""
class Single_Taxa(object):
	def __init__(self):
		""" 
		this variable signifies the cluster index that the taxa belongs 
		if the value is -1 (< 0) then it is still not inserted in a cluster 
		otherwise it is part of a valid cluster 
		"""
		self.clust_idx_part = -1

	def _Get_Taxa_Part_Clust_Idx(self):
		return self.clust_idx_part
		
	def _Set_Clust_Idx_taxa_Part(self, inp_clust_idx):
		self.clust_idx_part = inp_clust_idx
	
	# this function is called after formation of consensus tree
	def _PrintFinalTaxaInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n taxa key: ' + str(key))
		fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
		fp.close()

#-----------------------------------------------------
""" 
this class defines a couplet, according to the information obtained from input trees
key of this class --- taxa1, taxa2  
the consensus relation, frequencies, priority, and support scores 
in the class, the edge type signifies the relation between a pair of taxa
"""
class Reln_TaxaPair(object):
	def __init__(self):    
		""" 
		frequencies of individual relations 
		there are 4 types of edges (relationship) between a pair of taxa 
		"""
		self.freq_count = [0] * 4    
		""" 
		a connection priority value is defined as the 
		no of occurrences of this particular relation between this pair of taxa 
		minus the sum of no of occurrences of other relation types between this couplet
		"""
		self.priority_reln = [0] * 4    
		""" 
		this is the support score for different types of relations between a couplet
		"""
		self.support_score = [0] * 4
		""" 
		For this couplet, it stores the extra gene count with respect to all the gene trees
		"""
		self.XL_sum_gene_trees = []
		""" 
		this variable stores the no of trees supporting the taxa pair 
		"""
		self.supporting_trees = 0
		"""
		this list contains the union of taxa list underlying the LCA of this couplet
		for individual input trees
		"""
		self.LCA_Underlying_Taxa_List = []
		"""
		for a couplet xy, and for a relation r3, following array has 3 elements:
		1) count when level of x < level of y (relation r1 similar)
		2) count when level of x > level of y (relation r2 similar)
		3) count when level of x = level of y (relation r3 similar)
		"""
		self.ALL_Reln_Level_Diff_Info_Count = [0] * 3
		self.ALL_Reln_Level_Diff_Val_Count = [0] * 3
		
		"""
		this is a variable containing the binned average of the XL values
		of very high frequency
		initially the value is set as -1, to signify that the computation is not done
		once the computation (for a couplet) is done, the value is subsequently used and returned
		"""
		self.binned_avg_XL = -1
		
		"""
		this is a list containing the number of instances
		when R4 relation is actually a pseudo R1 relation
		"""
		self.freq_R4_pseudo_R1R2 = [0] * 2
		"""
		this is a list of allowed relations between this couplet
		we allow only those relations which have a significant frequency 
		compared to the number of supporting trees
		"""
		self.allowed_reln_list = []
		
		"""
		this list contains the sorted frequency count in descending order
		starting from the highest to lowest
		"""
		self.sorted_freq_count = [0] * 4

	#----------------------------------
	"""
	this function creates the sorted freq list in descending order
	starting from the highest to the lowest
	"""
	def _SortFreqList(self):
		for r in range(4):
			self.sorted_freq_count[r] = self.freq_count[r]
		# sort the list
		self.sorted_freq_count.sort(reverse=True)

	"""
	this function checks whether the input frequency value 
	is either max or 2nd max frequency
	from all different relations and corresponding frequency values
	"""
	def _CheckMaxor2ndMaxReln(self, inp_freq):
		"""
		as the sorted frequency list is stored in descending order
		so we check whether the given frequency lies in the 1st or 2nd index
		of course, frequency should be > 0
		"""
		if (inp_freq > 0):
			if (inp_freq == self.sorted_freq_count[0]) or (inp_freq == self.sorted_freq_count[1]):
				return 1
		return 0
		
	#----------------------------------
	"""
	this function returns the ratio of level value
	@param: idx: 	if 0, returns ratio w.r.t r1
							  if 1, returns ratio w.r.t r2
	"""
	def _GetLevelValRatio(self, idx):
		level_val_r1 = self.ALL_Reln_Level_Diff_Val_Count[0]
		level_val_r2 = self.ALL_Reln_Level_Diff_Val_Count[1]
		if (idx == 0):
			if ((level_val_r1 + level_val_r2) > 0):
				return (level_val_r1 * 1.0) / (level_val_r1 + level_val_r2)
			else:
				return 0
		else:	#if (idx == 1):
			if ((level_val_r1 + level_val_r2) > 0):
				return (level_val_r2 * 1.0) / (level_val_r1 + level_val_r2)
			else:
				return 0
			
		return 0
	#----------------------------------
	"""
	Appends underlying taxon set of the LCA node for this couplet
	"""
	def _AppendUnderlyingTaxonList(self, inp_list):
		self.LCA_Underlying_Taxa_List = list(set(self.LCA_Underlying_Taxa_List) | set(inp_list))

	"""
	Get the complete set of taxon, underlying the LCA nodes, for all input trees, 
	corresponding to this couplet
	"""
	def _GetUnderlyingTaxonList(self):
		return self.LCA_Underlying_Taxa_List

	#----------------------------------
	"""
	this function checks whether the 'target_reln' (RELATION_R1 or RELATION_R2) 
	can be established between this couplet
	"""
	def CheckHigherPriority(self, target_reln):
		if (target_reln == RELATION_R1):
			if (self._CheckTargetRelnConsensus(RELATION_R1) == True):
				if (self.ALL_Reln_Level_Diff_Info_Count[0] == max(self.ALL_Reln_Level_Diff_Info_Count)):
					if (self.ALL_Reln_Level_Diff_Val_Count[0] > self.ALL_Reln_Level_Diff_Val_Count[1]):
						return 1

		if (target_reln == RELATION_R2):
			if (self._CheckTargetRelnConsensus(RELATION_R2) == True):
				if (self.ALL_Reln_Level_Diff_Info_Count[1] == max(self.ALL_Reln_Level_Diff_Info_Count)):
					if (self.ALL_Reln_Level_Diff_Val_Count[1] > self.ALL_Reln_Level_Diff_Val_Count[0]):
						return 1
	
		return 0

	#----------------------------------
	"""
	these functions adjust the set of allowed relations among a couplet
	"""
	def _AddAllowedReln(self, inp_reln):
		self.allowed_reln_list.append(inp_reln)
		
	def _GetAllowedRelnList(self):
		return self.allowed_reln_list
	
	def _RemoveAllowedReln(self, inp_reln):
		if inp_reln in self.allowed_reln_list:
			self.allowed_reln_list.remove(inp_reln)
			
	#----------------------------------
	def _AddFreqPseudoR1(self, idx, r=1):
		self.freq_R4_pseudo_R1R2[idx] = self.freq_R4_pseudo_R1R2[idx] + r
		
	def _GetFreqPseudoR1(self, idx):
		return self.freq_R4_pseudo_R1R2[idx]
		
	def _NormalizeR1R2LevelDiff(self):
		for i in range(3):
			self.ALL_Reln_Level_Diff_Val_Count[i] = (self.ALL_Reln_Level_Diff_Val_Count[i] * 1.0) / self.supporting_trees
		
	#----------------------------------
	"""
	this function computes the level difference (relation based) count 
	# parameters:
	@src_reln: corresponds to the level of relation from which the subtraction will take place
	@dest_reln: if between 0 and 2, corresponds to the subtracted relation
	if it is -1, all relations except the src_reln will be subtracted
	@abs_comp: if true, will compute the absolute value of the subtracted quantity
	@norm: if true, will normalize the subtracted value with the number of supporting trees
	"""
	def _GetRelnLevelDiff(self, src_reln, dest_reln, abs_comp, norm):
		reln_array = [RELATION_R1, RELATION_R2, RELATION_R3]
		src_reln_idx = reln_array.index(src_reln)

		target_val = self.ALL_Reln_Level_Diff_Info_Count[src_reln_idx]
		
		if (dest_reln == -1):
			# all relations will be subtracted
			for dest_reln_idx in range(3):
				if (dest_reln_idx == src_reln_idx):
					continue
				target_val = target_val - self.ALL_Reln_Level_Diff_Info_Count[dest_reln_idx]
		else:
			dest_reln_idx = reln_array.index(dest_reln)
			target_val = target_val - self.ALL_Reln_Level_Diff_Info_Count[dest_reln_idx]
			
		if (abs_comp == True):
			target_val = math.fabs(target_val)
			
		if (norm == True):
			target_val = (target_val * 1.0) / self.supporting_trees
			
		return target_val

	"""
	this function checks whether the input relation type is a consensus relation
	among this couplet
	"""
	def _CheckTargetRelnConsensus(self, inp_reln_type):
		if (self.freq_count[inp_reln_type] == max(self.freq_count)):
			return True
		return False
	
	def _CheckTargetRelnLevelConsensus(self, src_reln):
		reln_array = [RELATION_R1, RELATION_R2, RELATION_R3]
		src_reln_idx = reln_array.index(src_reln)
		sum_level_count = sum(self.ALL_Reln_Level_Diff_Info_Count)
		sum_level_val = self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1]
		if (src_reln_idx == 0) or (src_reln_idx == 1):
			if (self.ALL_Reln_Level_Diff_Info_Count[src_reln_idx] > (0.5 * sum_level_count)) and \
				(self.ALL_Reln_Level_Diff_Val_Count[src_reln_idx] > (0.5 * sum_level_val)):
				return 1
		else:
			if (self.ALL_Reln_Level_Diff_Info_Count[2] > (self.ALL_Reln_Level_Diff_Info_Count[0] \
				+ self.ALL_Reln_Level_Diff_Info_Count[1])):
				return 1
			
		return 0
		
	def _IncrAllRelnLevelDiffInfoCount(self, idx, val):
		self.ALL_Reln_Level_Diff_Info_Count[idx] = self.ALL_Reln_Level_Diff_Info_Count[idx] + 1
		self.ALL_Reln_Level_Diff_Val_Count[idx] = self.ALL_Reln_Level_Diff_Val_Count[idx] + val
			
	#----------------------------------
	"""
	this function adds one supporting tree for this couplet
	"""
	def _AddSupportingTree(self):
		self.supporting_trees = self.supporting_trees + 1

	"""
	this function returns the number of input trees supporting this couplet
	"""
	def _GetNoSupportTrees(self):
		return self.supporting_trees
		
	#----------------------------------
	"""
	this function adds one XL value, computed for a particular input tree
	to the list of XL values for this couplet
	"""
	def _AddXLVal(self, val):
		self.XL_sum_gene_trees.append(val)
	
	def _GetXLList(self):
		return self.XL_sum_gene_trees
	
	def _GetXLSumGeneTrees(self):
		return sum(self.XL_sum_gene_trees)
		
	"""
	this function computes the average of XL measures
	"""
	def _GetAvgXLGeneTrees(self):
		return (sum(self.XL_sum_gene_trees) * 1.0) / self.supporting_trees

	"""
	function to return the average of XL values for this couplet
	depending on the user parameters, average, median, or binned average XL is returned
	"""
	def _GetNormalizedXLSumGeneTrees(self, dist_type):
		if (dist_type == 1):
			return self._GetAvgXLGeneTrees()
		elif (dist_type == 2):
			# average of mean and mode
			return (self._GetAvgXLGeneTrees() + self._GetMultiModeXLVal()) / 2.0

	#------------------------------------------------
	"""
	this function computes the binned average of XL values associated for this couplet
	"""
	def _GetMultiModeXLVal(self, Output_Text_File=None):
		if (self.binned_avg_XL == -1):
			
			Bin_Width = (1.0 / MODE_BIN_COUNT)
			len_list = [0] * MODE_BIN_COUNT
			
			if Output_Text_File is not None:
				fp = open(Output_Text_File, 'a') 
			
			# sort the XL list
			self.XL_sum_gene_trees.sort()
			
			for j in range(len(self.XL_sum_gene_trees)):
				curr_xl_val = self.XL_sum_gene_trees[j]
				bin_idx = int(curr_xl_val / Bin_Width)
				if (bin_idx == MODE_BIN_COUNT):
					bin_idx = bin_idx - 1
				len_list[bin_idx] = len_list[bin_idx] + 1
			
			if Output_Text_File is not None:
				for i in range(MODE_BIN_COUNT):
					fp.write('\n bin idx: ' + str(i) + ' len:  ' + str(len_list[i]))
			
			# this is the maximum length of a particular bin
			# corresponding to max frequency
			max_freq = max(len_list)
			
			if Output_Text_File is not None:
				fp.write('\n Max freq: ' + str(max_freq))
			
			num = 0
			denom = 0
			for i in range(MODE_BIN_COUNT):
				if (len_list[i] >= (MODE_PERCENT * max_freq)):
					list_start_idx = sum(len_list[:i])
					list_end_idx = list_start_idx + len_list[i] - 1
					value_sum = sum(self.XL_sum_gene_trees[list_start_idx:(list_end_idx+1)])
					num = num + value_sum
					denom = denom + len_list[i]
					if Output_Text_File is not None:
						fp.write('\n Included bin idx: ' + str(i) + ' starting point: ' + str(list_start_idx) \
							+ 'ending point: ' + str(list_end_idx) + ' sum: ' + str(value_sum))
			
			self.binned_avg_XL = (num / denom)
			
			if Output_Text_File is not None:
				fp.write('\n Final binned average XL: ' + str(self.binned_avg_XL))
				fp.close()
			
		return self.binned_avg_XL

	#----------------------------------
	"""
	this function returns the frequency of the consensus relation
	"""
	def _GetConsensusFreq(self):
		return max(self.freq_count)

	"""
	this function returns the frequency of the input relation
	specified by the variable 'reln_type'
	"""
	def _GetEdgeWeight(self, reln_type):
		return self.freq_count[reln_type]      
	
	"""
	this function adds a specified frequency count (default 1)
	corresponding to the relation specified by the variable 'reln_type'
	"""
	def _AddEdgeCount(self, reln_type, val=1):
		self.freq_count[reln_type] = self.freq_count[reln_type] + val
	
	#----------------------------------
	""" 
	this function computes the support score value associated with individual couplet
	for all different relations
	"""
	def _SetCostMetric(self):
		for reln_type in range(4):
			# assign the score metric for this relation type
			self.support_score[reln_type] = self.freq_count[reln_type] * self.priority_reln[reln_type]
	
	"""
	this function returns the support score of the input relation
	specified by the variable 'reln_type'
	"""
	def _GetEdgeCost_ConnReln(self, reln_type):
		return self.support_score[reln_type]
	
	"""
	this function updates (increments) the support score of the input relation
	specified by the variable 'reln_type'
	and by the amount 'incr_cost'
	"""
	def _IncrEdgeCost_ConnReln(self, reln_type, incr_cost):
		self.support_score[reln_type] = self.support_score[reln_type] + incr_cost
	
	#----------------------------------
	"""
	this function returns the priority value for a given input relation 
	specified by the variable 'reln_type'
	"""
	def _GetConnPrVal(self, reln_type):
		return self.priority_reln[reln_type]
	
	""" 
	this function calculates connection priority value for 
	each of the relation types
	"""
	def _SetConnPrVal(self):
		"""
		this is the sum of frequencies for all the relation types
		"""
		listsum = sum(self.freq_count)
		"""
		now determine the connection priority of a 
		particular relation type with respect to other relations     
		"""
		for reln_type in range(4):
			"""
			here we use the difference of current relation type 
			frequency with the frequencies of all other relations
			"""
			self.priority_reln[reln_type] = 2 * self.freq_count[reln_type] - listsum
		
	"""
	this function checks whether a couplet is non-conflicting
	that is, only one relation between them exists throughout all the gene trees
	in such a case, a binary variable 1 and the corresponding 
	relation type is returned in the form of a list
	otherwise, a value 0 and a defaukt relation R4 is returned
	"""
	def _CheckNonConflictingCouplet(self):
		# this is the sum of frequencies for all the relation types
		listsum = sum(self.freq_count)
		""" 
		this code section is used when there exists an unique relation (non-conflicting couplet)
		and we try to detect it
		"""
		outlist = [0, RELATION_R4]
		for reln_type in range(4):
			if (self.freq_count[reln_type] == listsum) and (listsum > 0):
				outlist = [1, reln_type]
				break
			elif (self.freq_count[reln_type] > 0) and (self.freq_count[reln_type] < listsum):
				break
		return outlist

	"""
	this function prints all the information associated with a couplet
	"""
	def _PrintRelnInfo(self, key, Output_Text_File):
		level_val_diff_R1 = self.ALL_Reln_Level_Diff_Val_Count[0]
		level_val_diff_R2 = self.ALL_Reln_Level_Diff_Val_Count[1]
		level_val_denom = (level_val_diff_R1 + level_val_diff_R2)
		level_val_diff = (level_val_diff_R1 - level_val_diff_R2)
		
		fp = open(Output_Text_File, 'a')    
		fp.write('\n\n\n taxa pair key: ' + str(key))
		fp.write('\n relations [type/count/priority_reln/score]: ')
		for i in range(4):
			fp.write('\n [' + str(i) + '/' \
				+ str(self.freq_count[i]) + '/' + str(self.priority_reln[i]) \
				+ '/' + str(self.support_score[i]) + ']')
		fp.write('\n Sorted freq list: ' + str(self.sorted_freq_count))
		fp.write('\n AVERAGE Sum of excess gene **** : ' + str(self._GetAvgXLGeneTrees()))
		#fp.write('\n Mode Sum of extra lineage **** : ' + str(self._GetMultiModeXLVal()))
		#fp.write('\n Mean(Avg + Mode) of extra lineage **** : ' \
			#+ str((self._GetAvgXLGeneTrees() + self._GetMultiModeXLVal()) / 2.0))
		fp.write('\n No of supporting trees : ' + str(self.supporting_trees))
		fp.write('\n ALL relation based Level diff info count (r1/r2/r3): ' \
			+ str(self.ALL_Reln_Level_Diff_Info_Count))
		fp.write('\n ALL relation based Level diff Val count (r1/r2/r3): ' \
			+ str(self.ALL_Reln_Level_Diff_Val_Count))
		fp.write('\n R4 relation pseudo (R1/R2) count: ' + str(self.freq_R4_pseudo_R1R2))
		fp.write('\n Allowed relation list: ' + str(self.allowed_reln_list))
		if (level_val_denom > 0):
			fp.write('\n Level Val r1 reln ratio : ' + str((level_val_diff_R1 * 1.0) / level_val_denom))
			fp.write('\n Level Val r2 reln ratio : ' + str((level_val_diff_R2 * 1.0) / level_val_denom))
			fp.write('\n Level Val r3/r4 reln ratio : ' + str((math.fabs(level_val_diff)) / level_val_denom))
		fp.close()

#-----------------------------------------------------
""" 
this class is representative of a cluster of taxa that are related via equality relationship 
according to the rule of equivalence partition 
"""
class Cluster_node(object):
	def __init__(self, inp_taxa=None):
		"""
		taxa list of the current cluster
		"""
		self.Species_List = [] 
		"""
		# set to 1 once the cluster is traversed during DFS order of traversing the clusters
		this is required in printing the supertree in newick format 
		"""
		self.explored = 0    
		"""
		stores the indices of clusters cy such that curr_clust->cy is achieved
		"""
		self.out_edge_list = []
		"""
		stores the indices of clusters cy such that cy->curr_clust is achieved
		"""
		self.in_edge_list = []
		"""
		stores the indices of clusters cy such that cy and curr_clust 
		are not related by any directed edge based connection
		"""
		self.no_edge_list = []
		"""
		stores the indices of clusters cy such that curr_clust->cy connection 
		needs to be checked
		"""
		self.possible_R1_list = []
		"""
		during initialization, append one tuple to this cluster
		"""
		if inp_taxa is not None:
			self._Append_taxa(inp_taxa)    

	#---------------------------------------------
	"""
	these functions keep track whether the cluster node is used 
	during newick string formation for supertree construction
	each of the clusters (containing a set of taxa) should be visited 
	exactly once for supertree generation
	"""
	def _SetExploredStatus(self):
		self.explored = 1

	def _ResetExploredStatus(self):
		self.explored = 0
		
	def _GetExploredStatus(self):
		return self.explored
	#---------------------------------------------
	"""
	following function modifies and returns (if required) 
	the constituent taxa list within this cluster
	"""
	def _GetSpeciesList(self):
		return self.Species_List
	
	def _GetCardinality(self):
		return len(self.Species_List)
				
	# append one species information in this cluster
	def _Append_taxa(self, inp_taxa):
		if inp_taxa not in self.Species_List:
			self.Species_List.append(inp_taxa)
	#---------------------------------------------
	"""
	returns the number of clusters cy such that cy->curr_clust connection is present
	"""
	def _Get_Indegree(self):
		return len(self.in_edge_list)

	"""
	returns the number of clusters cy such that curr_clust->cy connection is present
	"""
	def _Get_Outdegree(self):
		return len(self.out_edge_list)
	
	"""
	returns the list of clusters cy such that curr_clust->cy connection is present
	"""
	def _GetOutEdgeList(self):
		return self.out_edge_list
	
	"""
	returns the list of clusters cy such that cy->curr_clust connection is present
	"""
	def _GetInEdgeList(self):
		return self.in_edge_list    

	"""
	returns the set of clusters cy such that there exists no directed edge 
	between cy and curr_clust is present
	"""
	def _GetNoEdgeList(self):
		return self.no_edge_list    
	
	"""
	adds one cluster cy in the list of clusters such that curr_clust->cy connection is present
	"""
	def _AddOutEdge(self, dest_clust_idx):
		if dest_clust_idx not in self.out_edge_list:
			self.out_edge_list.append(dest_clust_idx)
		
	"""
	adds one cluster cy in the list of clusters such that cy->curr_clust connection is present
	"""
	def _AddInEdge(self, src_clust_idx):
		if src_clust_idx not in self.in_edge_list:
			self.in_edge_list.append(src_clust_idx)

	"""
	adds one cluster cy in the list of clusters such that there exists no directed edge 
	between cy and curr_clust
	"""
	def _AddNoEdge(self, src_clust_idx):
		if src_clust_idx not in self.no_edge_list:
			self.no_edge_list.append(src_clust_idx)
			
	"""
	removes one cluster cy in the list of clusters such that curr_clust->cy connection is now deleted
	"""
	def _RemoveOutEdge(self, dest_clust_idx):
		if dest_clust_idx in self.out_edge_list:
			self.out_edge_list.remove(dest_clust_idx)    
		
	"""
	removes one cluster cy in the list of clusters such that cy->curr_clust connection is now deleted
	"""
	def _RemoveInEdge(self, dest_clust_idx):
		if dest_clust_idx in self.in_edge_list:
			self.in_edge_list.remove(dest_clust_idx)    
		
	"""
	removes one cluster cy in the list of clusters such that information between no 
	directed edge connection between the curr_clust and the cluster cy is now deleted
	"""
	def _RemoveNoEdge(self, dest_clust_idx):
		if dest_clust_idx in self.no_edge_list:
			self.no_edge_list.remove(dest_clust_idx)    
			
	#--------------------------------------------------------
	# add - sourya
	def _AddPossibleR1(self, dest_clust_idx):
		if dest_clust_idx not in self.possible_R1_list:
			self.possible_R1_list.append(dest_clust_idx)
	
	def _RemovePossibleR1(self, dest_clust_idx):
		if dest_clust_idx in self.possible_R1_list:
			self.possible_R1_list.remove(dest_clust_idx)
			
	def _GetPossibleR1List(self):
		return self.possible_R1_list
	#--------------------------------------------------------
		
	def _PrintClusterInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n cluster key: ' + str(key))
		fp.write('\n species list: ' + str(self.Species_List))
		#print 'its indegree: ', self.indegree
		#print 'its outdegree: ', self.outdegree
		fp.write('\n out edge list: ' + str(self.out_edge_list))
		fp.write('\n in edge list: ' + str(self.in_edge_list))
		fp.write('\n Possible R1 list: ' + str(self.possible_R1_list))
		fp.close()    
