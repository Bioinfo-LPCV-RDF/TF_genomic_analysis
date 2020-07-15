# -*- coding: utf-8 -*-

### ----------------------------
### TFBS Interdistances program
### ----------------------------


'''
This program allows to calculate interdistances between transcription factor binding sites.
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks), and unbound sequences).
You can display a plot for several thresholds.
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import re
import time
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from interdistances_functions import *
from plot_functions import *
from operator import truediv
import operator
from collections import Counter
import matplotlib.patches as mpatches
from matplotlib import pylab
import types
import argparse
import logging
from optparse import OptionParser
from scipy import stats
import pickle
parser = argparse.ArgumentParser()                                               


parser.add_argument("--positive", "-pos", nargs='*')
parser.add_argument("--negative_sets", "-neg", nargs='*')
parser.add_argument("--legend", "-leg", nargs='*')
parser.add_argument("--matrix_type", "-mt", type=str, default="freq")
parser.add_argument("--threshold", "-th",nargs='+',type = float, default= [-8, -9, -10])
parser.add_argument("--Interdistance_maxValue", "-maxInter", type = int, default= 20)
parser.add_argument("--Interdistance_minValue", "-minInter", type = int, default= 0)
parser.add_argument("--histo", "-histo", type = bool, default= False)
parser.add_argument("--points", "-points", type = bool, default= True)
parser.add_argument("--curve", "-curve", type = bool, default= False)
parser.add_argument("--sum_threshold", "-sum_threshold", type = bool, default= False)
parser.add_argument("--save", "-save",action='store_true', default= False)
parser.add_argument("--load", "-load",action='store_true', default= False)
parser.add_argument("--write_inter", "-write_inter",action='store_true', default= False)
parser.add_argument("--no_threshold_panel", "-no_th",action='store_true', default= False)
parser.add_argument("--one_panel", "-one_panel",action='store_true', default= False)
parser.add_argument("--no_absolute_panel", "-no_absolute_panel",action='store_true', default= False)
parser.add_argument("--log", "-log",action='store_true', default= False)
parser.add_argument("--matrix", "-mat",type=str, default="")
parser.add_argument("--pseudoCount", "-pc",type = float, default = 0.0)
parser.add_argument("--dependancies", "-d", type=str, default="")
parser.add_argument("--tffm", "-tffm", type=str, default="")
parser.add_argument("--output", "-o", type=str, default="./")
parser.add_argument("--name", "-n", type=str, default="figure.svg")
parser.add_argument("--offset_left", "-ol", type = int, default=0)
parser.add_argument("--offset_right", "-or", type = int, default=0)
parser.add_argument("--max_y_axis", "-maxy", type = float, default=0)

args = parser.parse_args()

negative_sets = args.negative_sets
print(negative_sets)
legend = args.legend
save = args.save
load = args.load
write_inter = args.write_inter
threshold = args.threshold
Interdistance_maxValue = args.Interdistance_maxValue
Interdistance_minValue = args.Interdistance_minValue
#Interdistance_minValue=0 ## this feature doees not work yet
histo = args.histo
points = args.points
sum_threshold = args.sum_threshold
positive=args.positive
one_panel=args.one_panel
log=args.log
no_threshold_panel=args.no_threshold_panel
no_absolute_panel=args.no_absolute_panel
maxy=args.max_y_axis


if len(positive)!=1 and load is False:
	print("list of positive sets and load = False non compatible, exiting")
	quit()
if not load:
	positive=positive[0]

tffm=args.tffm
matrix=args.matrix
#either matrix ot tffm
if tffm == "":
        matrixArg=args.matrix
        matrixType=args.matrix_type
elif tffm != "" and matrix != "":
        print("Chose either tffm or matrix (pfm), exiting")
	quit()
elif tffm != "":
        import sys
        sys.path.append("/home/304.3-STRUCTPDEV/lib/TFFM")
        import tffm_module
        from constants import TFFM_KIND
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

if (Interdistance_minValue != 0 and load is True) or (Interdistance_minValue != 0 and save is True) :
        print("load,save options only work when minInter=0. Exiting")
	quit()

        


pseudoCount = args.pseudoCount
dependencyFile = args.dependancies
figure_name = args.name
output_directory = args.output
offset_right = args.offset_right
offset_left = args.offset_left
output= output_directory +"/"+ figure_name

neg_flag=False # can be true if save option called
pos_flag=False # can be true if save option called

if histo == True and len(threshold) > 1 : 
	print("Impossible to display an histogram with several thresholds, you have to change a parameter.")
	sys.exit(0)
	
########################################### About the main matrix (true if tffm="")  #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas...

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''

################################# If save == True  #################################
if save and load:
	print("Save = load = True, non compatible options, exiting")
	quit()

if save :
	basename_pos_file=os.path.basename(os.path.splitext(positive)[0])
	basename_neg_file=os.path.basename(os.path.splitext(negative_sets[0])[0])
	list_th=''
	for elt in threshold:
		list_th=list_th+'_'+str(elt)
	if os.path.exists(output_directory+'/'+basename_neg_file+list_th+'_neg.pkl'):
		neg_flag=True
	if os.path.exists(output_directory+'/'+basename_pos_file+list_th+'_pos.pkl'):
		pos_flag=True
	if pos_flag and neg_flag:
		print("Save = True, The required files already exist, exiting")
		quit()
################################# If load == True  #################################

if load:
	norm=list()
	if len(positive)==1:
		comp=False
		positive=positive[0]
		with open(positive,'rb') as f1:
			pos_pickler=pickle.Unpickler(f1)
			scores_pos=pos_pickler.load()
			InterDR = scores_pos['DR']
			InterER = scores_pos['ER']
			InterIR = scores_pos['IR']
			threshold_pos = scores_pos['threshold']
			maxInter_pos = scores_pos['maxInter']
			len_pos = scores_pos['len_pos']
		with open(negative_sets[0],'rb') as f1:
			neg_pickler=pickle.Unpickler(f1)
			scores_neg=neg_pickler.load()
			InterDR_N = scores_neg['DR_N']
			InterER_N = scores_neg['ER_N']
			InterIR_N = scores_neg['IR_N']
			threshold_neg = scores_neg['threshold']
			maxInter_neg = scores_neg['maxInter']
			len_neg = scores_neg['len_neg']
		norm=[len_pos/len_neg]*len(threshold)
		if threshold_pos != threshold_neg or maxInter_neg != maxInter_pos:
			print("load=True, Files not corresponding to the same parameters, exiting")
			quit()
	else:
		comp=True
		list_pos=list()
		InterDR = []
		InterER = []
		InterIR = []
		InterDR_N = []
		InterER_N = []
		InterIR_N = []
		for th in threshold:
			for elt in positive:
				with open(elt,'rb') as f1:
					pos_pickler=pickle.Unpickler(f1)
					scores_pos=(pos_pickler.load())
					threshold_pos=scores_pos['threshold']
					which=threshold_pos.index(th)
					InterDR.append(scores_pos['DR'][which])
					InterER.append(scores_pos['ER'][which])
					InterIR.append(scores_pos['IR'][which])
					maxInter_pos = scores_pos['maxInter']
					len_pos = scores_pos['len_pos']
				with open(negative_sets[0],'rb') as f1:
					neg_pickler=pickle.Unpickler(f1)
					scores_neg=(neg_pickler.load())
					threshold_neg=scores_neg['threshold']
					which=threshold_neg.index(th)
					InterDR_N.append(scores_neg['DR_N'][which])
					InterER_N.append(scores_neg['ER_N'][which])
					InterIR_N.append(scores_neg['IR_N'][which])
					maxInter_neg = scores_neg['maxInter']
					len_neg = scores_neg['len_neg']
				norm.append(len_pos/len_neg)
				if threshold_pos != threshold_neg or maxInter_neg != maxInter_pos:
					print("load=True, Files not corresponding to the same parameters, exiting")
					quit()
		threshold=legend

	Interdistance_maxValue=maxInter_pos
	interdist_sum = []
	interdist_sum_N = []
	rate = []	

	for a,b,c,d,e,f,g in zip(InterDR,InterER,InterIR,InterDR_N,InterER_N,InterIR_N,norm) :
		interdist_sum.append(sum(a) + sum(b) + sum(c))
		interdist_sum_N.append(sum(d) + sum(e) + sum(f))
		rate.append([divide(sum(a), sum(d))/g, divide(sum(b) , sum(e))/g, divide(sum(c) , sum(f))/g])


	relative_DR = []
	relative_ER = []
	relative_IR = []
	relative_DR_neg = []
	relative_ER_neg = []
	relative_IR_neg = []

	for a,b,c,d,e in zip(threshold,InterDR,interdist_sum,InterDR_N,interdist_sum_N) :
		relative_DR.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_DR_neg.append( [divide(x , float(e)) for x in d] )
	for a,b,c,d,e in zip(threshold,InterER,interdist_sum,InterER_N,interdist_sum_N) :
		relative_ER.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_ER_neg.append( [divide(x , float(e)) for x in d] )
	for a,b,c,d,e in zip(threshold,InterIR,interdist_sum,InterIR_N,interdist_sum_N) :
		relative_IR.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_IR_neg.append( [divide(x , float(e)) for x in d] )
	command = sys.argv			
        print comp
        if points == True and negative_sets and one_panel == False :
	        negative_sets_points(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command,output,comp,log)
        if points == True and negative_sets and one_panel == True :
	        negative_sets_points_ER(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command,output,comp,log,no_threshold_panel)

	quit()


####################################################################################
FastaFile=positive
if tffm =="":
        MatrixFile=matrixArg
        # These 3 lines allows to retrieve the matrix from the file
        F = open(MatrixFile,"r")
        matrix = F.read().replace("\r","\n")
        print matrix
        F.close()
        # These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list
        num = re.compile(r"([+-]?\d+[.,]\d+)")
        Mdata = num.findall(matrix)
        #print("list(reversed(Mdata)) : ",list(reversed(Mdata)))
        matScore, lenMotif = get_score_matrix(Mdata,matrixType,pseudoCount)
        # The following line allows to produce the reversed matrix
        '''if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
        Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
        So we can calculate with this reverse matrix, the score of the complementary strand.
        '''
        matRev = list(reversed(matScore))
        tffm_first_order =""
else:
        matRev=""
        matScore=""
        tffm_first_order = tffm_module.tffm_from_xml(tffm, TFFM_KIND.FIRST_ORDER)
        lenMotif=tffm_first_order.__len__() + 2
	
########## get INTERDISTANCE VALUES for POSITIVE sets:
if not pos_flag:
	len_pos=0
	sequence_number=0

	with open(FastaFile,"r") as f1:
		for line in f1:
			if line.find(">") != -1:	
				line=line.strip()
				line=line.replace("-",":")
				line=line.split(":")
				len_pos+=float(line[2]) - float(line[1]) + 1 - lenMotif
				sequence_number+=1
	len_pos=float(len_pos)
	print(sequence_number)
	sequence_number_pos=sequence_number

	InterDR, InterER, InterIR = get_interdist(matScore,matRev,FastaFile,threshold,offset_left,offset_right,Interdistance_maxValue,sum_threshold,lenMotif,dependencyFile,sequence_number,tffm_first_order,write_inter,output)
	##### Create empty lists to store interdistances occurences for the negative set:

	InterDR_N = []
	InterER_N = []
	InterIR_N = []
	lenThr = 0
	listThr = []
	for a in threshold :
		InterDR_N.append( [0] * (Interdistance_maxValue + 1) )
		InterER_N.append( [0] * (Interdistance_maxValue + 1) )
		InterIR_N.append( [0] * (Interdistance_maxValue + 1) )	
		listThr.append(lenThr)
		lenThr = lenThr + 1


########## get INTERDISTANCE occurences for NEGATIVE sets			       
if not neg_flag:
	len_neg=0
	if negative_sets :
		i=0
		for fastafileN in negative_sets :
			i+=1
			sequence_number=0
			with open(fastafileN,"r") as f1:
				for line in f1:
					if line.find(">") != -1:
						line=line.strip()
						line=line.replace("-",":")
						line=line.split(":")
						len_neg+=float(line[2])-float(line[1]) + 1 - lenMotif
						sequence_number+=1
			len_neg=float(len_neg)
			print(sequence_number)

			sequence_number_neg=sequence_number              

			InterDR_N_temp, InterER_N_temp, InterIR_N_temp = get_interdist(matScore,matRev,fastafileN,threshold,offset_left,offset_right,Interdistance_maxValue,sum_threshold,lenMotif,dependencyFile,sequence_number,tffm_first_order,write_inter,output+"_neg"+str(i))
			# addition of the occurences of every negative sets
			for a,b,c,d in zip(InterDR_N_temp,InterER_N_temp,InterIR_N_temp,listThr) :
				InterDR_N[d] = [x + y for x, y in zip(InterDR_N[d], a)]
				InterER_N[d] = [x + y for x, y in zip(InterER_N[d], b)]
				InterIR_N[d] = [x + y for x, y in zip(InterIR_N[d], c)]
		# divide by the number of negative sets in order to have a mean
		if len(negative_sets) > 0 :
			for a,b,c,d in zip(InterDR_N,InterER_N,InterIR_N,listThr) :
				InterDR_N[d] = [x / float(1) for x in a]
				InterER_N[d] = [x / float(1) for x in b]
				InterIR_N[d] = [x / float(1) for x in c]


############## this block does not work yet, for now, Interdistance_minValue is set to 0 ###########
# remove the distances from 0 to Interdistance_minValue
ind=0
while ind < len(threshold):
	InterDR[ind]=InterDR[ind][Interdistance_minValue: len(InterDR[ind])]
	InterER[ind]=InterER[ind][Interdistance_minValue: len(InterER[ind])]
	InterIR[ind]=InterIR[ind][Interdistance_minValue: len(InterIR[ind])]
	if len(negative_sets) > 0 :                                
		InterDR_N[ind]=InterDR_N[ind][Interdistance_minValue: len(InterDR_N[ind])]
		InterER_N[ind]=InterER_N[ind][Interdistance_minValue: len(InterER_N[ind])]
		InterIR_N[ind]=InterIR_N[ind][Interdistance_minValue: len(InterIR_N[ind])]
	ind+=1
####################################################################################################        
        
if save :
	if not os.path.exists(output_directory+'/'+basename_pos_file+list_th+'_pos.pkl'):
		with open(output_directory+'/'+basename_pos_file+list_th+'_pos.pkl','wb') as f1:
			scores_pos={
				"len_pos": len_pos,
				"DR": InterDR,
				"ER": InterER,
				"IR": InterIR,
				"threshold": threshold,
				"maxInter": Interdistance_maxValue,
				}
			pickler_neg=pickle.Pickler(f1)
			pickler_neg.dump(scores_pos)
	if not os.path.exists(output_directory+'/'+basename_neg_file+list_th+'_neg.pkl'):
		with open(output_directory+'/'+basename_neg_file+list_th+'_neg.pkl','wb') as f1:
			scores_neg={
				"len_neg": len_neg,
				"DR_N": InterDR_N,
				"ER_N": InterER_N,
				"IR_N": InterIR_N,
				"threshold": threshold,
				"maxInter": Interdistance_maxValue,
				}
			pickler_neg=pickle.Pickler(f1)
			pickler_neg.dump(scores_neg)
	quit()



norm=len_pos/len_neg

print "norm = " + str(norm)

interdist_sum = []
interdist_sum_N = []
rate = []


if negative_sets :
	for a,b,c,d,e,f in zip(InterDR,InterER,InterIR,InterDR_N,InterER_N,InterIR_N) :
		interdist_sum.append(sum(a) + sum(b) + sum(c))
		interdist_sum_N.append(sum(d) + sum(e) + sum(f))
		rate.append([divide(sum(a), sum(d))/norm, divide(sum(b) , sum(e))/norm, divide(sum(c) , sum(f))/norm])
		
print ("rate")
print rate


relative_DR = []
relative_ER = []
relative_IR = []
relative_DR_neg = []
relative_ER_neg = []
relative_IR_neg = []

k=0
if negative_sets :	
	k+=1
	for a,b,c,d,e in zip(threshold,InterDR,interdist_sum,InterDR_N,interdist_sum_N) :
		relative_DR.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_DR_neg.append( [divide(x , float(e)) for x in d] )
	for a,b,c,d,e in zip(threshold,InterER,interdist_sum,InterER_N,interdist_sum_N) :
		relative_ER.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_ER_neg.append( [divide(x , float(e)) for x in d] )
	for a,b,c,d,e in zip(threshold,InterIR,interdist_sum,InterIR_N,interdist_sum_N) :
		relative_IR.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_IR_neg.append( [divide(x , float(e)) for x in d] )
	if one_panel:
		relative_ER = []
		relative_ER_neg = []
		for a,b,c,d,e,f,g,h,i in zip(threshold,InterER,interdist_sum,InterER_N,interdist_sum_N, InterDR, InterDR_N, InterIR, InterIR_N) :
			relative_ER.append( [divide(x+y+z , float(c)) for x,y,z in zip(b,h,f)] )
			if len(negative_sets) > 0 :
				relative_ER_neg.append( [divide(x+y+z , float(e)) for x,y,z in zip(d,i,g)] )
#print("relative_DR : ",relative_DR)
#print("relative_ER : ",relative_ER)
#print("relative_IR : ",relative_IR)
if write_inter:
	maxi=[]
	i=0
	if one_panel:
		for a,b in zip(relative_ER,relative_ER_neg):
			maxi.append([])
			enrichment=map(divide, a, b)
			maxi[i].append(enrichment)
			maxi[i].append(enrichment)
			maxi[i].append(enrichment)
			i+=1
	else:
		for a,b,c,d,e,f in zip(relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg):
			maxi.append([])
			maxi[i].append(map(divide, a, b))
			maxi[i].append(map(divide, c, d))
			maxi[i].append(map(divide, e, f))
			i+=1
	with open(output_directory+'/enrichments',"w") as OUT:
		i=0
		for elt in maxi:
			OUT.write(str(threshold[i])+"\n")
			i+=1
			for elt2 in elt:
				line=""
				for elt3 in elt2:
					if line !="":
						line+="\t"+str(elt3)
					else:
						line=str(elt3)
				OUT.write(line+"\n")

#if write_inter:
#	maxi=[]
#	c=0
#	if one_panel:
#		while c<=2:
#			maxi.append([])
#			for a,b in zip(relative_ER,relative_ER_neg) :
#				maxi[c].append(map(divide, a, b))
#			c+=1
#	else:
#		maxi.append([])
#		for a,b in zip(relative_DR,relative_DR_neg) :
#			maxi[c].append(map(divide, a, b))
#		c+=1
#		maxi.append([])
#		for a,b in zip(relative_ER,relative_ER_neg) :
#			maxi[c].append(map(divide, a, b))
#		c+=1
#		maxi.append([])
#		for a,b in zip(relative_IR,relative_IR_neg) :
#			maxi.append(map(divide, a, b))
#		
#	with open(output_directory+'/enrichments',"w") as OUT:
#		for elt in maxi:
#			for elt2 in elt:
#				line=""
#				for elt3 in elt2:
#					if line !="":
#						line+="\t"+str(elt3)
#					else:
#						line=str(elt3)
#				OUT.write(line+"\n")

command = sys.argv			
			
			
if histo == True and negative_sets :
	negative_sets_histo(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command)
	
if histo == False and points == False and negative_sets :
	negative_sets_curve(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command)

if points == True and negative_sets and one_panel == False :
	negative_sets_points(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command,output,Interdistance_minValue,False,log)

if points == True and negative_sets and one_panel == True :
	output= output_directory +"/"+ figure_name
	negative_sets_points_ER(Interdistance_maxValue,relative_ER,relative_ER_neg,threshold,rate,command,output,Interdistance_minValue,False,log,no_threshold_panel,no_absolute_panel,maxy)
	#output= output_directory +"/IR_"+ figure_name
	#negative_sets_points_ER(Interdistance_maxValue,relative_IR,relative_IR_neg,threshold,rate,command,output,Interdistance_minValue,False,log,no_threshold_panel,no_absolute_panel,maxy)
	#output= output_directory +"/DR_"+ figure_name
	#negative_sets_points_ER(Interdistance_maxValue,relative_DR,relative_DR_neg,threshold,rate,command,output,Interdistance_minValue,False,log,no_threshold_panel,no_absolute_panel,maxy)
