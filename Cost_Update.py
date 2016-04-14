#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

##-----------------------------------------
#"""
#this function checks whether the couplet (t1,t2) is supported and have 'target_reln' as its consensus 
#or whether the couplet (t2,t1) is supported and have 'Complementary_Reln(target_reln)' as its consensus 
#@Return: 
#success / failure: 1 / 0
#couplet_key: key of the couplet
#target_reln: the consensus relation between this couplet
#"""
#def Check_Couplet_With_Consensus_Reln(t1, t2, target_reln):
	#key1 = (t1, t2)
	#key2 = (t2, t1)
	#if (key1 in TaxaPair_Reln_Dict):
		#couplet_key = key1
		#if (TaxaPair_Reln_Dict[couplet_key]._CheckTargetRelnConsensus(target_reln) == True):
			#return 1, couplet_key, target_reln
	#elif (key2 in TaxaPair_Reln_Dict):
		#couplet_key = key2
		#if (TaxaPair_Reln_Dict[couplet_key]._CheckTargetRelnConsensus(Complementary_Reln(target_reln)) == True):
			#return 1, couplet_key, Complementary_Reln(target_reln)
	
	#return 0, None, RELATION_R4

##-----------------------------------------
#"""
#this function checks a set of three consensus relations 
#which are mutually conflicting
#the couplets are (A,B), (B,C), and (C,A)
#in such a case, it removes one of the consensus relations
#"""
#def Check_Consensus_Reln_Conflict():
	#for t1 in COMPLETE_INPUT_TAXA_LIST:
		#for t2 in COMPLETE_INPUT_TAXA_LIST:
			#if (t1 == t2):
				#continue
			#"""
			#check whether the couplet (t1, t2) has a consensus relation of R3
			#or if the couplet (t2, t1) has a consensus relation of R3
			#"""
			#c1_support_cons, c1_key, c1_reln = Check_Couplet_With_Consensus_Reln(t1, t2, RELATION_R3)
			#if (c1_support_cons == 0):
				#continue
			#c1_support_score = TaxaPair_Reln_Dict[c1_key]._GetEdgeCost_ConnReln(c1_reln)
			#c1_freq = TaxaPair_Reln_Dict[c1_key]._GetEdgeWeight(c1_reln)
			
			#for t3 in COMPLETE_INPUT_TAXA_LIST:
				#if (t1 == t3) or (t2 == t3):
					#continue
				
				#"""
				#check whether the couplet (t2, t3) has a consensus relation of R1
				#or if the couplet (t3, t2) has a consensus relation of R2
				#"""
				#c2_support_cons, c2_key, c2_reln = Check_Couplet_With_Consensus_Reln(t2, t3, RELATION_R1)
				#if (c2_support_cons == 0):
					#continue
				#c2_support_score = TaxaPair_Reln_Dict[c2_key]._GetEdgeCost_ConnReln(c2_reln)
				#c2_freq = TaxaPair_Reln_Dict[c2_key]._GetEdgeWeight(c2_reln)
				
				#"""
				#check whether the couplet (t3, t1) has a consensus relation of R3
				#or if the couplet (t1, t3) has a consensus relation of R3
				#"""
				#c3_support_cons, c3_key, c3_reln = Check_Couplet_With_Consensus_Reln(t3, t1, RELATION_R3)
				#if (c3_support_cons == 0):
					#continue
				#c3_support_score = TaxaPair_Reln_Dict[c3_key]._GetEdgeCost_ConnReln(c3_reln)
				#c3_freq = TaxaPair_Reln_Dict[c3_key]._GetEdgeWeight(c3_reln)
				
				#"""
				#here we have obtained the desired couplets
				#which are basically partitions of a triplet
				#along with their conflicting consensus relations
				#here we have to remove the relation "c1_reln" and the corresponding couplet "c1_key"
				#from the list of support scores in the priority queue
				#"""
				#if (c1_freq > c3_freq):
					#target_support_score = c1_support_score
					#target_sublist = [c1_key[0], c1_key[1], c1_reln, target_support_score]
					#if target_sublist in Cost_List_Taxa_Pair_Multi_Reln:
						#Cost_List_Taxa_Pair_Multi_Reln.remove(target_sublist)
					#if target_sublist in Cost_List_Taxa_Pair_Single_Reln:
						#Cost_List_Taxa_Pair_Single_Reln.remove(target_sublist)
  
#------------------------------------------------------
""" 
this code section implements the max priority queue 
"""
#------------------------------------------------------
# parent node of current node
def Parent(idx):
	return int((idx-1) / 2)

# left child node of current node  
def Left(idx):
	return (2*idx+1)

# right child node of current node
def Right(idx):
	return (2*idx+2)

#------------------------------------------------------------------------
# version 3 - latest - sourya - 6th April 2016
#-------------------------
"""
this function is used for cost based sorting of the inp_queue in the unweighted supertree
this function compares two elements of the heap
returns 1 if element in the i'th index is lower than element in the j'th index
that is, j'th index has higher priority than i'th index
"""
def Lower_Score_Value(inp_queue, i, j):
	key1 = (inp_queue[i][0], inp_queue[i][1])
	reln1 = inp_queue[i][2]
	score1 = inp_queue[i][3]

	key2 = (inp_queue[j][0], inp_queue[j][1])
	reln2 = inp_queue[j][2]
	score2 = inp_queue[j][3]

	if (score1 > 0) or (score2 > 0):
		if (score1 < score2):
			return 1
		elif (score1 > score2):
			return 0
		elif (score1 == score2):
			"""
			following checks are performed in tie case of cost
			"""
			if (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) < TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
				return 1      
			elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) > TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
				return 0
			elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == True) \
				and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == False):
				return 0
			elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == False) \
				and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == True):
				return 1
			#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) < TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
				#return 1
			#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) > TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
				#return 0    
			#elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
				#return 0
			#elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
				#return 1
			#elif (reln1 == RELATION_R3) and (reln2 != RELATION_R3):
				#return 1
			#elif (reln1 != RELATION_R3) and (reln2 == RELATION_R3):
				#return 0
			#elif (((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 == RELATION_R3) or (reln2 == RELATION_R4))):
				#return 0
			#elif (((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 == RELATION_R3) or (reln1 == RELATION_R4))):
				#return 1
			#elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
				#return 0
			#elif ((reln2 == RELATION_R4) and (reln1 == RELATION_R3)):
				#return 1
	else:
		"""
		here both scores are either zero or negative
		so we directly check the frequency and priority measures
		without comparing the score values
		"""
		if (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) < TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
			return 1      
		elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) > TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
			return 0
		elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == True) \
			and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == False):
			return 0
		elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == False) \
			and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == True):
			return 1
		#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) < TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
			#return 1
		#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) > TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
			#return 0    
		#elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
			#return 0
		#elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
			#return 1
		#elif (reln1 == RELATION_R3) and (reln2 != RELATION_R3):
			#return 1
		#elif (reln1 != RELATION_R3) and (reln2 == RELATION_R3):
			#return 0
		#elif (((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 == RELATION_R3) or (reln2 == RELATION_R4))):
			#return 0
		#elif (((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 == RELATION_R3) or (reln1 == RELATION_R4))):
			#return 1
		#elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
			#return 0
		#elif ((reln2 == RELATION_R4) and (reln1 == RELATION_R3)):
			#return 1

	return 0

##------------------------------------------------------------------------
## version 2 - this is the latest - sourya - 1st nov 2015
##-------------------------
#"""
#this function is used for cost based sorting of the inp_queue in the unweighted supertree
#this function compares two elements of the heap
#returns 1 if element in the i'th index is lower than element in the j'th index
#that is, j'th index has higher priority than i'th index
#"""
#def Lower_Score_Value(inp_queue, i, j):
	#key1 = (inp_queue[i][0], inp_queue[i][1])
	#reln1 = inp_queue[i][2]
	#score1 = inp_queue[i][3]

	#key2 = (inp_queue[j][0], inp_queue[j][1])
	#reln2 = inp_queue[j][2]
	#score2 = inp_queue[j][3]

	#if (score1 > 0) or (score2 > 0):
		#if (score1 < score2):
			#return 1
		#elif (score1 > score2):
			#return 0
		#elif (score1 == score2):
			#"""
			#following checks are performed in tie case of cost
			#"""
			#if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == True) \
				#and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == False):
				#return 0
			#elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == False) \
				#and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == True):
				#return 1
			#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) < TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
				#return 1
			#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) > TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
				#return 0    
			#elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) < TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
				#return 1      
			#elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) > TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
				#return 0
			#elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
				#return 0
			#elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
				#return 1
			#elif (reln1 == RELATION_R3) and (reln2 != RELATION_R3):
				#return 1
			#elif (reln1 != RELATION_R3) and (reln2 == RELATION_R3):
				#return 0
			##elif (((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 == RELATION_R3) or (reln2 == RELATION_R4))):
				##return 0
			##elif (((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 == RELATION_R3) or (reln1 == RELATION_R4))):
				##return 1
			##elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
				##return 0
			##elif ((reln2 == RELATION_R4) and (reln1 == RELATION_R3)):
				##return 1
	#else:
		#"""
		#here both scores are either zero or negative
		#so we directly check the frequency and priority measures
		#without comparing the score values
		#"""
		#if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == True) \
			#and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == False):
			#return 0
		#elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(reln1) == False) \
			#and (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(reln2) == True):
			#return 1
		#elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) < TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
			#return 1      
		#elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) > TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
			#return 0
		#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) < TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
			#return 1
		#elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) > TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
			#return 0    
		#elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
			#return 0
		#elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
			#return 1
		#elif (reln1 == RELATION_R3) and (reln2 != RELATION_R3):
			#return 1
		#elif (reln1 != RELATION_R3) and (reln2 == RELATION_R3):
			#return 0
		##elif (((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 == RELATION_R3) or (reln2 == RELATION_R4))):
			##return 0
		##elif (((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 == RELATION_R3) or (reln1 == RELATION_R4))):
			##return 1
		##elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
			##return 0
		##elif ((reln2 == RELATION_R4) and (reln1 == RELATION_R3)):
			##return 1

	#return 0

#--------------------------------------------------------------------
# this function exchanges two elements in the heap
def Exchange_Elem(inp_queue, i, j):
	temp_key = inp_queue[i][0]
	inp_queue[i][0] = inp_queue[j][0]
	inp_queue[j][0] = temp_key
	temp_key = inp_queue[i][1]
	inp_queue[i][1] = inp_queue[j][1]
	inp_queue[j][1] = temp_key
	temp_edge = inp_queue[i][2]
	inp_queue[i][2] = inp_queue[j][2]
	inp_queue[j][2] = temp_edge
	temp_val = inp_queue[i][3]
	inp_queue[i][3] = inp_queue[j][3]
	inp_queue[j][3] = temp_val

# maintain max heap property
# note: heap_size may not be the actual length of the queue
# but the working length (on which the remaining sorting operation will take place)
def Max_Heapify(inp_queue, idx, heap_size):
	l = Left(idx)
	r = Right(idx)
	if (l < heap_size) and Lower_Score_Value(inp_queue, idx, l):
		largest_idx = l
	else:
		largest_idx = idx
	if (r < heap_size) and Lower_Score_Value(inp_queue, largest_idx, r):
		largest_idx = r
	if (largest_idx != idx):
		Exchange_Elem(inp_queue, idx, largest_idx)
		Max_Heapify(inp_queue, largest_idx, heap_size)

# extract the current maximum and also pop the element from the heap
def Heap_Extract_Max(inp_queue):
	if (len(inp_queue) < 1):
		print 'underflow of max priority queue'
	# 1st element is the maximum
	max_elem = list(inp_queue[0])
	# replace the first element with the last element of the queue
	Exchange_Elem(inp_queue, 0, len(inp_queue) - 1)
	# delete the last element of the queue
	del inp_queue[len(inp_queue) - 1]
	heap_size = len(inp_queue)
	# call the max_heapify function to maintain the heap property
	# 0 is the starting index of the list storing the heap structure
	Max_Heapify(inp_queue, 0, heap_size)
	return max_elem
  
# this function builds the priority queue (max heap property)
def Build_Max_Heap(inp_queue):
	heap_size = len(inp_queue)
	for idx in range(int(len(inp_queue) / 2), -1, -1):
		Max_Heapify(inp_queue, idx, heap_size)
    
# this is the heap sort algorithm
def Sort_Priority_Queue(inp_queue):
	Build_Max_Heap(inp_queue)
	heap_size = len(inp_queue)
  
