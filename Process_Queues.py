#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import Conflict_Detect
from Conflict_Detect import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------------------
""" 
this function processes support score queue containing couplet based relations
"""
def Proc_Queue_Couplet(Reachability_Graph_Mat, Output_Text_File):
	Inp_Queue = Queue_Score_Couplet
	
	while (0 < len(Inp_Queue)):
		""" 
		extract the 1st element of "Inp_Queue" 
		since it is sorted to have max cost at the beginning 
		"""
		outlist = Heap_Extract_Max(Inp_Queue)
		
		src_taxa_idx = outlist[0]
		src_taxa_label = COMPLETE_INPUT_TAXA_LIST[src_taxa_idx]
		dest_taxa_idx = outlist[1]
		dest_taxa_label = COMPLETE_INPUT_TAXA_LIST[dest_taxa_idx]
		reln_type = outlist[2]
		reln_freq = outlist[3]
		conn_score = outlist[4]
		
		"""
		there is no concept of conflict in this case
		we just check whether R3 relation is predominant among all taxa pairs within this pair of cluster
		"""
		clust1 = Taxa_Info_Dict[src_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		clust2 = Taxa_Info_Dict[dest_taxa_idx]._Get_Taxa_Part_Clust_Idx()
		
		if (DEBUG_LEVEL > 0):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ==>>>>>>>>> current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
					' Clust1: ' + str(clust1) + ' Clust2: ' + str(clust2) + \
						' relation type: ' + str(reln_type) + ' relation freq: ' + str(reln_freq) + \
							' conn score: ' + str(conn_score))
			fp.close()

		conflict_detection = Possible_Conflict_Curr_Reln(clust1, clust2, Reachability_Graph_Mat, reln_type, Output_Text_File)
		if (conflict_detection == 0):
			if (DEBUG_LEVEL > 0):
				fp = open(Output_Text_File, 'a')    
				fp.write('\n ==>>>>>>>>> NEW CONN ---- relation type: ' + str(reln_type) \
					+ '  from ' + str(src_taxa_label) + ' to ' + str(dest_taxa_label) + \
						' - clust: ' + str(clust1) + ' to clust: ' + str(clust2) \
						+ ' relation freq: ' + str(reln_freq) + ' conn score: ' + str(conn_score))
				fp.close()
			"""
			also update the reachability graph information
			"""
			Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, clust1, clust2, reln_type, Output_Text_File)

	return Reachability_Graph_Mat

