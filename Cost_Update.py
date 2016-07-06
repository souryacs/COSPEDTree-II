#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#----------------------------------------------------------------
"""
here we fill the support score queues with the couplet relations and the support score values
from which the DAG will be constructed
"""
def Fill_Support_Score_Queues_Couplet_Based():
	""" 
	here, we process all the couplets
	individual couplets lead to individual cluster pairs
	for individual relations between a couplet, corresponding edge between the cluster pair is established
	"""  
	for l in TaxaPair_Reln_Dict:
		#r1_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R1)
		#r2_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R2)
		#r3_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R3)
		#r4_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(RELATION_R4)

		#"""
		#single_edge_exist: if TRUE, means that only one type of relation 
		#is supported (with respect to input trees) between this couplet
		#detection of it during setting the priority values of different relations
		#basically, we are looking for the consensus relation
		#"""
		#single_edge_exist_list = TaxaPair_Reln_Dict[l]._CheckNonConflictingCouplet()
		#single_edge_exist = single_edge_exist_list[0]
		#consensus_reln_type = single_edge_exist_list[1]
		
		"""
		maximum frequency among all 4 relations
		"""
		max_freq = TaxaPair_Reln_Dict[l]._GetConsensusFreq()

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

		for reln_type in [RELATION_R4, RELATION_R1, RELATION_R2, RELATION_R3]:
			curr_reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type)
			if (curr_reln_freq > 0):
				TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
				sublist = [l[0], l[1], reln_type, TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type), \
					TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(reln_type)]
				Queue_Score_Couplet.append(sublist)
				"""
				also add the consensus relation in a separate list
				"""
				if (curr_reln_freq == max_freq):
					TaxaPair_Reln_Dict[l]._AddConsensusRelnList(reln_type)
				
		#"""
		#for other types of couplets, not included in above criterion
		#"""
		#for reln_type in [RELATION_R4, RELATION_R1, RELATION_R2, RELATION_R3]:
			#curr_reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(reln_type)
			#if (curr_reln_freq > 0):
				#if (reln_type == RELATION_R4):
					#"""
					#non zero frequency of R4 relation is considered
					#"""
					#TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
				#elif (reln_type == RELATION_R1) or (reln_type == RELATION_R2):
					#if (r4_freq == 0) or (curr_reln_freq >= (0.7 * max_freq)):	# or (curr_reln_freq == max_freq_nonR3):
					##if (curr_reln_freq >= (0.7 * max_freq)) or (curr_reln_freq == max_freq_nonR3):
						#"""
						#either zero frequency of R4 relation, or sufficient frequency of R1/R2 relations are considered
						#for their inclusion in the 
						#"""
						#TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
				#else:	#if (reln_type == RELATION_R3):
					#if (len(TaxaPair_Reln_Dict[l]._GetAllowedRelnList()) <= 1) \
						#and (TaxaPair_Reln_Dict[l]._CheckTargetRelnConsensus(reln_type) == True):
						#"""
						#first add relation R3 in the allowed relation list
						#"""
						#TaxaPair_Reln_Dict[l]._AddAllowedReln(reln_type)
						#"""
						#differentiate between two cases
						#1) When R3 is the sole relation between this couplet
						#2) When R3 is the majority consensus relation between this couplet
						#"""
						#if (consensus_reln_type == RELATION_R3) and (single_edge_exist == 1):
							#"""
							#we also add R1, R2, R4 relations as the allowed relations
							#"""
							#TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R1)
							#TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R2)
							#TaxaPair_Reln_Dict[l]._AddAllowedReln(RELATION_R4)
							#"""
							#add the support score relation in the queue
							#"""
							#reln_freq = TaxaPair_Reln_Dict[l]._GetEdgeWeight(consensus_reln_type)
							##support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(consensus_reln_type)
							#sublist = [l[0], l[1], consensus_reln_type, reln_freq, reln_freq]	#, support_score]
							#Queue_Score_R3_SingleReln.append(sublist)
						#elif (TaxaPair_Reln_Dict[l]._CheckTargetRelnMajorityConsensus(RELATION_R3) == True) and (r3_freq > 1):
							#"""
							#we impose one condition that the no of trees bearing this R3 relation should be > 1
							#add the support score relation in the queue
							#"""
							##support_score = TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(RELATION_R3)
							#sublist = [l[0], l[1], RELATION_R3, r3_freq, r3_freq]	#, support_score]
							#Queue_Score_R3_MajCons.append(sublist)
	
	
	return

#------------------------------------------------------
""" this code section implements the max priority queue """
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

##------------------------------------------------------------------------
## version 1X - 15th Sep, 2015 
##-------------------------
#"""
#this function is used for sorting the support score queue (priority queue structure)

#Parameters:
#1) inp_queue: Input priority queue
	#can be either Nonconflicting support score queue or the conflicting support score queue
#2) i, j denotes indices of two different elements within the priority queue
#these two elements are to be compared

#output:
#returns index with the higher support score measure
#"""
#def Higher_Score_Value(inp_queue, i, j):
	#key1 = (inp_queue[i][0], inp_queue[i][1])
	#reln1 = inp_queue[i][2]
	#freq1 = inp_queue[i][3]
	#score1 = inp_queue[i][4]

	#key2 = (inp_queue[j][0], inp_queue[j][1])
	#reln2 = inp_queue[j][2]
	#freq2 = inp_queue[j][3]
	#score2 = inp_queue[j][4]

	#"""
	#function when non-conflicting support score queue is used
	#"""
	## case A - if both scores are strictly positive
	#if (score1 > 0) or (score2 > 0):
		#if (score1 < score2):
			#return j
		#elif (score1 > score2):
			#return i
		#else:	#if (score1 == score2):
			#if (freq1 < freq2):
				#return j      
			#elif (freq2 < freq1):
				#return i
			##elif ((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 != RELATION_R1) and (reln2 != RELATION_R2)):
				##return i
			##elif ((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 != RELATION_R1) and (reln1 != RELATION_R2)):
				##return j
			##elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
				##return i
			##elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
				##return j
	#else:
		#"""
		#here one or both scores are either zero or negative
		#"""
		#if (freq1 < freq2):
			#return j      
		#elif (freq2 < freq1):
			#return i     
		##elif ((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 != RELATION_R1) and (reln2 != RELATION_R2)):
			##return i
		##elif ((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 != RELATION_R1) and (reln1 != RELATION_R2)):
			##return j
		##elif (reln1 == RELATION_R4) and (reln2 != RELATION_R4):
			##return i
		##elif (reln1 != RELATION_R4) and (reln2 == RELATION_R4):
			##return j
	
	## default return condition
	#return i
	
#------------------------------------------------------------------------

def Resolve_R4_and_R4(inp_queue, i, j):
	key_i = (inp_queue[i][0], inp_queue[i][1])
	reln_i = inp_queue[i][2]
	freq_i = inp_queue[i][3]
	score_i = inp_queue[i][4]
	priority_i = TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i)
	
	key_j = (inp_queue[j][0], inp_queue[j][1])
	reln_j = inp_queue[j][2]
	freq_j = inp_queue[j][3]
	score_j = inp_queue[j][4]
	priority_j = TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)

	# first we prioritize the consensus relation
	if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return i
	elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return j
	# next we prioritize the relation with higher frequency
	elif (freq_i < freq_j):
		return j
	elif (freq_i > freq_j):
		return i
	# next we prioritize the relation with higher priority
	elif (priority_i < priority_j):
		return j
	elif (priority_i > priority_j):
		return i  
	# last, we select the relation with high sum of XL
	elif (TaxaPair_Reln_Dict[key_i]._GetAvgXLGeneTrees() < TaxaPair_Reln_Dict[key_j]._GetAvgXLGeneTrees()):
		return j
	else:
		return i

def Resolve_R4_and_R3(inp_queue, i, j):
	key_i = (inp_queue[i][0], inp_queue[i][1])
	reln_i = inp_queue[i][2]
	freq_i = inp_queue[i][3]
	score_i = inp_queue[i][4]
	priority_i = TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i)
	
	key_j = (inp_queue[j][0], inp_queue[j][1])
	reln_j = inp_queue[j][2]
	freq_j = inp_queue[j][3]
	score_j = inp_queue[j][4]
	priority_j = TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)

	# first we prioritize the consensus relation
	if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return i
	elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return j
	# next we prioritize the relation with higher frequency
	elif (freq_i < freq_j):
		return j
	elif (freq_i > freq_j):
		return i 
	# next we prioritize the relation with higher priority
	elif (priority_i < priority_j):
		return j
	elif (priority_i > priority_j):
		return i  
	else:	#if (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
		# same couplets - prioritize R4
		return i
	#else:
		## for different couplet, prioritize R3
		#return j

def Resolve_R4_and_R1_R2(inp_queue, i, j):
	key_i = (inp_queue[i][0], inp_queue[i][1])
	reln_i = inp_queue[i][2]
	freq_i = inp_queue[i][3]
	score_i = inp_queue[i][4]
	priority_i = TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i)
	
	key_j = (inp_queue[j][0], inp_queue[j][1])
	reln_j = inp_queue[j][2]
	freq_j = inp_queue[j][3]
	score_j = inp_queue[j][4]
	priority_j = TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)

	# first we prioritize the consensus relation
	if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return i
	elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return j
	# next we prioritize the relation with higher frequency  
	elif (freq_i < freq_j):
		return j
	elif (freq_i > freq_j):
		return i
	# next we prioritize the relation with higher priority
	elif (priority_i < priority_j):
		return j
	elif (priority_i > priority_j):
		return i  
	elif (inp_queue[i][0] == inp_queue[j][0]) and (inp_queue[i][1] == inp_queue[j][1]):
		# same couplets - prioritize R4
		return i
	else:
		## for different couplets, prioritize R1 / R2 relation
		#return j
		# last, we select the relation with high sum of XL
		if (TaxaPair_Reln_Dict[key_i]._GetAvgXLGeneTrees() < TaxaPair_Reln_Dict[key_j]._GetAvgXLGeneTrees()):
			return j
		else:
			return i
  
def Resolve_R3_and_R1_R2(inp_queue, i, j):
	key_i = (inp_queue[i][0], inp_queue[i][1])
	reln_i = inp_queue[i][2]
	freq_i = inp_queue[i][3]
	score_i = inp_queue[i][4]
	priority_i = TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i)
	
	key_j = (inp_queue[j][0], inp_queue[j][1])
	reln_j = inp_queue[j][2]
	freq_j = inp_queue[j][3]
	score_j = inp_queue[j][4]
	priority_j = TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)

	# first we prioritize the consensus relation
	if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return i
	elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return j
	# next we prioritize the relation with higher frequency    
	elif (freq_i < freq_j):
		return j
	elif (freq_i > freq_j):
		return i 
	# next we prioritize the relation with higher priority
	elif (priority_i < priority_j):
		return j
	elif (priority_i > priority_j):
		return i

	# R3 is set to low priority
	return j  

def Resolve_R3_and_R3(inp_queue, i, j):
	key_i = (inp_queue[i][0], inp_queue[i][1])
	reln_i = inp_queue[i][2]
	freq_i = inp_queue[i][3]
	score_i = inp_queue[i][4]
	priority_i = TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i)
	
	key_j = (inp_queue[j][0], inp_queue[j][1])
	reln_j = inp_queue[j][2]
	freq_j = inp_queue[j][3]
	score_j = inp_queue[j][4]
	priority_j = TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)

	# first we prioritize the consensus relation
	if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return i
	elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return j
	# next we prioritize the relation with higher frequency    
	elif (freq_i < freq_j):
		return j
	elif (freq_i > freq_j):
		return i 
	# next we prioritize the relation with higher priority
	elif (priority_i < priority_j):
		return j
	elif (priority_i > priority_j):
		return i    

	return i

def Resolve_R1_R2_and_R2_R1(inp_queue, i, j):
	key_i = (inp_queue[i][0], inp_queue[i][1])
	reln_i = inp_queue[i][2]
	freq_i = inp_queue[i][3]
	score_i = inp_queue[i][4]
	priority_i = TaxaPair_Reln_Dict[key_i]._GetConnPrVal(reln_i)
	
	key_j = (inp_queue[j][0], inp_queue[j][1])
	reln_j = inp_queue[j][2]
	freq_j = inp_queue[j][3]
	score_j = inp_queue[j][4]
	priority_j = TaxaPair_Reln_Dict[key_j]._GetConnPrVal(reln_j)

	# first we prioritize the consensus relation
	if (reln_i in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j not in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return i
	elif (reln_i not in TaxaPair_Reln_Dict[key_i]._GetConsensusRelnList()) and (reln_j in TaxaPair_Reln_Dict[key_j]._GetConsensusRelnList()):
		return j  
	# next we prioritize the relation with higher frequency    
	elif (freq_i < freq_j):
		return j
	elif (freq_i > freq_j):
		return i
	# next we prioritize the relation with higher priority
	elif (priority_i < priority_j):
		return j
	elif (priority_i > priority_j):
		return i    
	# last, we select the relation with high sum of XL
	elif (TaxaPair_Reln_Dict[key_i]._GetAvgXLGeneTrees() < TaxaPair_Reln_Dict[key_j]._GetAvgXLGeneTrees()):
		return j

	return i

#-------------------------------------------------

def Higher_Score_Value(inp_queue, i, j):
	key1 = (inp_queue[i][0], inp_queue[i][1])
	reln1 = inp_queue[i][2]
	freq1 = inp_queue[i][3]
	score1 = inp_queue[i][4]

	key2 = (inp_queue[j][0], inp_queue[j][1])
	reln2 = inp_queue[j][2]
	freq2 = inp_queue[j][3]
	score2 = inp_queue[j][4]
	
	if (score1 < score2):
		return j
	elif (score1 > score2):
		return i
	else:	#if (score1 == score2):
		# for tie case of cost
		# case A - if both relations are R4
		if ((reln1 == RELATION_R4) and (reln2 == RELATION_R4)):
			return Resolve_R4_and_R4(inp_queue, i, j)
		elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
			return Resolve_R4_and_R3(inp_queue, i, j)
		elif ((reln1 == RELATION_R3) and (reln2 == RELATION_R4)):
			return Resolve_R4_and_R3(inp_queue, j, i)
		elif ((reln1 == RELATION_R4) and ((reln2 == RELATION_R1) or (reln2 == RELATION_R2))):
			return Resolve_R4_and_R1_R2(inp_queue, i, j)
		elif ((reln2 == RELATION_R4) and ((reln1 == RELATION_R1) or (reln1 == RELATION_R2))):
			return Resolve_R4_and_R1_R2(inp_queue, j, i)
		elif ((reln1 == RELATION_R3) and ((reln2 == RELATION_R1) or (reln2 == RELATION_R2))):
			return Resolve_R3_and_R1_R2(inp_queue, i, j)
		elif ((reln2 == RELATION_R3) and ((reln1 == RELATION_R1) or (reln1 == RELATION_R2))):
			return Resolve_R3_and_R1_R2(inp_queue, j, i)
		elif (reln1 == RELATION_R3) and (reln2 == RELATION_R3):
			return Resolve_R3_and_R3(inp_queue, i, j)
		elif ((reln1 == RELATION_R1) and (reln2 == RELATION_R2)) or ((reln1 == RELATION_R2) and (reln2 == RELATION_R1)): 
			return Resolve_R1_R2_and_R2_R1(inp_queue, i, j)
		
	return i

#------------------------------------------------------------------------
"""
this function exchanges two elements in the heap
"""
def Exchange_Elem(inp_queue, i, j):
	"""
	here r corresponds to different fields (of the sublists)
	of the entries in the support score queue
	as individual queue entry consists of 5 different fields
	the copy operation should take care of the 5 fields
	"""
	total_elem = 5
	
	# swap individual fields
	for r in range(total_elem):
		temp_key = inp_queue[i][r]
		inp_queue[i][r] = inp_queue[j][r]
		inp_queue[j][r] = temp_key

	return

#-----------------------------------------------
"""
maintain max heap property
note: heap_size may not be the actual length of the queue
but the working length (on which the remaining sorting operation will take place)
"""
def Max_Heapify(inp_queue, idx, heap_size):
	l = Left(idx)
	r = Right(idx)
	if (l < heap_size) and (Higher_Score_Value(inp_queue, idx, l) == l):
		largest_idx = l
	else:
		largest_idx = idx
	if (r < heap_size) and (Higher_Score_Value(inp_queue, largest_idx, r) == r):
		largest_idx = r
	if (largest_idx != idx):
		Exchange_Elem(inp_queue, idx, largest_idx)
		Max_Heapify(inp_queue, largest_idx, heap_size)

#-----------------------------------------------
"""
extract the current maximum and also pop the element from the heap
"""
def Heap_Extract_Max(inp_queue):
	if (len(inp_queue) < 1):
		print 'underflow of max priority queue'
	"""
	1st element is the maximum
	"""
	max_elem = list(inp_queue[0])
	"""
	replace the first element with the last element of the queue
	"""
	Exchange_Elem(inp_queue, 0, len(inp_queue) - 1)
	"""
	delete the last element of the queue
	"""
	del inp_queue[len(inp_queue) - 1]
	heap_size = len(inp_queue)
	"""
	call the max_heapify function to maintain the heap property
	0 is the starting index of the list storing the heap structure
	"""
	Max_Heapify(inp_queue, 0, heap_size)
	return max_elem

#-----------------------------------------------
"""
this function builds the priority queue (max heap property)
"""
def Build_Max_Heap(inp_queue):
	heap_size = len(inp_queue)
	for idx in range(int(len(inp_queue) / 2), -1, -1):
		Max_Heapify(inp_queue, idx, heap_size)

#-----------------------------------------------
"""
this is the heap sort algorithm
parameters:
1) inp_queue: Input priority queue
	can be either Nonconflicting support score queue or the conflicting support score queue
"""
def Sort_Priority_Queue(inp_queue):
	Build_Max_Heap(inp_queue)
	heap_size = len(inp_queue)

