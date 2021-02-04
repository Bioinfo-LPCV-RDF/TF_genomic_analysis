# -*- coding: utf-8 -*-
#!usr/bin/python

### ----------------------------
### TFBS Interdistances program
### ----------------------------

######################				General intels				######################				
'''
This program allows to calculate interdistances between transcription factor binding sites.
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks), and unbound sequences).
You can display a plot for several thresholds.
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''
#python_version  :2.7.9
######################				import packages				######################				

import argparse
import numpy as np
import sys
import math
import re

######################				Functions				######################
def divide(a, b):
	try:
		return a/float(b)
	except:
		return 0.0

def compute_Length_seq(lib,len_motif):
	headers = np.unique(np.sort(lib['header']))
	length_set=0
	for header in headers:
		lengthseq=int(header.split(":")[-1].split("-")[1])-int(header.split(":")[-1].split("-")[0]) + 1 - len_motif
		length_set+=lengthseq
	return length_set

def define_spacing(lib, spacing_max, spacing_min, offset_l, offset_r, len_motif, conformations, thresholds, write_inter, output,suffix):
	
	headers = np.unique(np.sort(lib['header'])) # Getting all sequences to analyses.
	
	# preparing dictionnary to contain spacing results
	if len(conformations.keys()) == 0:
		for th in thresholds:
			conformations[th]={}
			conformations[th]['ER']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['IR']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['DR']=[0 for elt in range(0,spacing_max+1) ]
			conformations[th]['sum']=0
			conformations[th]['sumER']=0
			conformations[th]['sumIR']=0
			conformations[th]['sumDR']=0
	
	if write_inter:
		"""
		In case we need to report every spacing & conformation found.
		matricePosition1 & 2 are start position of the matrix for given score
		correctedPosition1 & 2 are start position of spacing distance (matrice position corrected by length of motif and possible offsets)
		"""
		with open(output+"_spacing_"+suffix+".tsv","w") as OUT:
			OUT.write("Peak" + "\t" + "Spacing" + "\t" + "Score1" + "\t" + "Score2" + "\t" + "matricePosition1" + "\t" + "matricePosition2" + "\t" + "correctedPosition1" + "\t" + "correctedPosition2" + "\n")
	minth=min(thresholds)
	for header in headers: # For each sequences
		sites = lib[np.where(lib['header'] == header)]
		sites = np.sort(sites,order='position')
		if len(sites) > 1: # if sequences contain at least 2 sites
			for x1 in range(0,len(sites)-1):
				elt1=sites[x1]
				for x2 in range(x1+1,len(sites)):
					elt2=sites[x2]
					"""
					double for loop to analyses possibles combinaisons of Spacing.
					DR: define by same strandness of both sites
					ER: first sites' strand must be forward and second sites' strand must be reverse
					IR: first sites' strand must be reverse and second sites' strand must be forward
					
					note that distance computation is dependant on conformation (offset used are not the same)
					
					"""
					conformation="NA"
					dist = "NA"
					if elt1["strand"] == elt2["strand"]: # DR
						if elt1["strand"] == "1" or elt1["strand"] == "+": #DR > >
							dist = (elt2["position"] + offset_l) - (elt1["position"] + len_motif - offset_r)
							position1=elt1["position"] + len_motif - offset_r
							position2=elt2["position"] + offset_l
						else: #DR < <
							dist = (elt2["position"] + offset_r) - (elt1["position"] + len_motif - offset_l)
							position1=elt1["position"] + len_motif - offset_l
							position2=elt2["position"] + offset_r
						conformation = "DR"
					if (elt1["strand"] == "-1" and elt2["strand"] =="1") or (elt1["strand"] == "-" and elt2["strand"] =="+"): #IR < >
						dist = (elt2["position"] + offset_r) - (elt1["position"] + len_motif - offset_r)
						conformation = "IR"
						position1=elt1["position"] + len_motif - offset_r
						position2=elt2["position"] + offset_r
					if (elt1["strand"] == "1" and elt2["strand"] =="-1") or (elt1["strand"] == "+" and elt2["strand"] =="-"): #ER > <
						dist = (elt2["position"] + offset_l) - (elt1["position"] + len_motif - offset_l)
						conformation = "ER"
						position1=elt1["position"] + len_motif - offset_l
						position2=elt2["position"] + offset_l
					if dist!= "NA" and (spacing_min <= dist <= spacing_max):
						for th in thresholds:
							if elt1['score'] > th and elt2['score'] > th: # making sure each sites have score higher than thresholds and add to corresponding sums
								conformations[th][conformation][dist]+=1
								conformations[th]['sum']+=1
								conformations[th]['sum'+str(conformation)]+=1
								if write_inter and th==minth: # report for every spacing/conformation found
									# report is done only for min threshold as you should be able to filter afterwards
									with open(output+"_spacing_"+suffix+".tsv","a") as OUT:
										OUT.write(header + "\t" + conformation + "_" + str(dist) + "\t" + str(elt1['score']) + "\t" + str(elt2['score']) + "\t" + str(elt1['position']) + "\t" + str(elt2['position']) + "\t" + str(position1) + "\t" + str(position2) + "\n")
	#for th in thresholds: # Debug report
		#print(str(th)+" - sum: "+str(conformations[th]['sum'])+"\nsum DR: "+str(conformations[th]['sumDR'])+"\nsumER: "+str(conformations[th]['sumER'])+"\nsumIR: "+str(conformations[th]['sumIR']))
	#for th in thresholds:
		#print("DR"+str(conformations[th]["DR"]))
		#print("ER"+str(conformations[th]["ER"]))
		#print("IR"+str(conformations[th]["IR"]))
	return conformations

def write_enrichment(relative_conformations, spacing_max, spacing_min, thresholds, output, no_absolute_panel, rate):
	dist_list = np.arange(spacing_min,spacing_max + 1)
	with open(output+"_enrichment.tsv","w") as OUT:
		OUT.write("threshold\tconformation\tspace\tenrichment\n")
		for th in thresholds:
			for conformation in ['DR', 'IR', 'ER']:
				for dist in dist_list:
					to_print=True
					if relative_conformations['neg'][th][conformation][dist]==0:
						print("======\nthreshold: "+str(th)+"\nconformation: "+str(conformation)+"\ndistance: "+str(dist)+"\nNo sites found in this conformation for negative set\nThreshold set too high\n=====")
						to_print=False
					if relative_conformations['pos'][th][conformation][dist]==0:
						print("======\nthreshold: "+str(th)+"\nconformation: "+str(conformation)+"\ndistance: "+str(dist)+"\nNo sites found in this conformation for positive set\nThreshold set too high\n=====")
						to_print=False
					if to_print:
						OUT.write(str(th) + "\t" + conformation + "\t" + str(dist) + "\t" + str(relative_conformations['pos'][th][conformation][dist]/relative_conformations['neg'][th][conformation][dist])+"\n")
			for dist in dist_list:
				OUT.write(str(th) + "\t" + 'all' + "\t" + str(dist) + "\t" + str(relative_conformations['pos'][th]['ALL'][dist]/relative_conformations['neg'][th]['ALL'][dist])+"\n")
	if no_absolute_panel:
		pass
	else:
		with open(output+"_rate.tsv","w") as OUT:
			OUT.write("threshold\tconformation\trate\n")
			i=0
			for th in thresholds:
				for conformation,indices in zip(['DR', 'IR', 'ER'],[0,2,1]):
					OUT.write(str(th) + "\t" + conformation + "\t" + str(rate[i][indices])+"\n")
				i+=1


######################				MAIN				######################
def main(	positive_file,
			negative_files,
			thresholds,
			spacing_max,
			spacing_min,
			offset_l,
			offset_r,
			len_motif,
			output,
			write_inter,
			no_absolute_panel
		):
	th_min = min(thresholds)
	#print(thresholds) # Debug print
	tffm=False
	# if at least one th is greater than 0, then tffm matrix have been used
	for th in thresholds:
		if th >0:
			tffm=True
	if tffm:
		positive_set = np.genfromtxt(positive_file, dtype=[('header','S300'), ('position',int), ('stop',int), ('strand','S5'), ('seq','S300'), ('type','S10'),('lengthM',int), ('score',float)])
		len_motif = (positive_set['stop'][1] - positive_set['position'][1] + 1)
		"""
		length of motif is extended for maximum number of sites computation because TFFM matrices needs 2 bp in front of the matrix to take in consideration probability of apparition of the site.
		We don't need this correction for spacing computation since score are reported with correct length of motif (length similar to PWM's length)
		"""
		length_positive=compute_Length_seq(positive_set, len_motif+2)
	else:
		positive_set = np.genfromtxt(positive_file, dtype=[('header','S300'),('position',int),('strand','S5'),('score',float)])
		length_positive=compute_Length_seq(positive_set, len_motif)
	
	# we filter out sequences & sites with score < minimum threhold to avoid unecessary computation (gain of time)
	#print("number of sequences to analyse before filtration by threshold: " + str(len(np.unique(np.sort(positive_set['header'])))))
	positive_set = positive_set[np.where( positive_set['score'] > th_min )]
	#print("number of sequences to analyse after filtration by threshold: " + str(len(np.unique(np.sort(positive_set['header'])))))
	conformations={}
	
	# finding ER, DR and IR conformation for each spacing
	print("defining spacing for positive set")
	conformations['pos']={}
	conformations['pos'] = define_spacing(positive_set, spacing_max, spacing_min, offset_l, offset_r, len_motif,conformations['pos'], thresholds, write_inter, output,"pos")
	
	
	conformations['neg']={}
	i=0
	length_negative=0
	for elt in negative_files:
		i+=1
		print("defining spacing for negative set "+str(i))
		if tffm:
			negative_set = np.genfromtxt(elt, dtype=[('header','S300'), ('position',int), ('stop',int), ('strand','S5'), ('seq','S300'), ('type','S10'),('lengthM',int), ('score',float)])
			"""
			length of motif is extended for maximum number of sites computation. See explanation done for positive set )
			"""
			length_negative+=compute_Length_seq(negative_set, len_motif+2)
		else:
			negative_set = np.genfromtxt(elt, dtype=[('header','S300'),('position',int),('strand','S5'),('score',float)])
			length_negative+=compute_Length_seq(negative_set, len_motif)
		negative_set = negative_set[np.where( negative_set['score'] > th_min )]
		#print("number of sequences to analyse after filtration by threshold: " + str(len(np.unique(np.sort(negative_set['header'])))))
		# finding ER, DR and IR conformation for each spacing & saving in one dictionnary for all negative set.
		conformations['neg'] = define_spacing(negative_set, spacing_max, spacing_min, offset_l, offset_r, len_motif, conformations['neg'], thresholds, write_inter, output,"neg"+str(i))
	# normalisation computation to take into account the higher number of sequences treated for the negative set
	norm = length_positive/float(length_negative)
	print("normalisation: "+str(norm))
	#print(conformations) # Debug print
	
	rate=[]
	for th in thresholds:
		"""
		Absolute Enrichment computation:
		( sum(DRpos) / sum(DRneg) ) * all_possible_site_pos/all_possible_site_neg
		"""
		rate.append(	[divide(conformations['pos'][th]['sumDR'],conformations['neg'][th]['sumDR'])/norm,
						divide(conformations['pos'][th]['sumER'],conformations['neg'][th]['sumER'])/norm,
						divide(conformations['pos'][th]['sumIR'],conformations['neg'][th]['sumIR'])/norm ])
	#print(rate) # Debug print
	
	relative_conformations={}
	relative_conformations['pos']={}
	relative_conformations['neg']={}
	"""
	Relative Enrichment computation: 
	sum(DRpos th|dist) / ( sum(NRpos) ) / sum(DRneg th|dist) / ( sum(NRneg) )
	
	NR : ER + DR + IR
	
	"""
	for key in relative_conformations.keys():
		for th in thresholds:
			relative_conformations[key][th]={}
			relative_conformations[key][th]['ER']=[ divide(x, float(conformations[key][th]['sum'])) for x in conformations[key][th]['ER'] ]
			relative_conformations[key][th]['IR']=[ divide(x, float(conformations[key][th]['sum'])) for x in conformations[key][th]['IR'] ]
			relative_conformations[key][th]['DR']=[ divide(x, float(conformations[key][th]['sum'])) for x in conformations[key][th]['DR'] ]
			relative_conformations[key][th]['ALL']=[ divide(x+y+z, float(conformations[key][th]['sum'])) for x,y,z in zip(conformations[key][th]['DR'],conformations[key][th]['ER'],conformations[key][th]['IR']) ]
	#for th in thresholds: # Debug prints
		# for conformation in ['DR', 'ER','IR','ALL']
			#print(conformation+" pos relative enrichment "+str(relative_conformations['pos'][th][conformation]))
			#print(conformation+" neg relative enrichment "+str(relative_conformations['neg'][th][conformation]))
	write_enrichment(relative_conformations, spacing_max, spacing_min, thresholds, output, no_absolute_panel, rate)
	
######################				Parsing arguments				######################

parser = argparse.ArgumentParser()
parser.add_argument("--positive", "-pos", help='scores file of seqences bound by TF, generated by scores.py')
parser.add_argument("--negative_sets", "-neg", nargs='*', help='list of scores files of seqences not bound by TF, generated by scores.py' )
parser.add_argument("--threshold", "-th",nargs='+',type = float, default= [-7, -8, -9, -10])
parser.add_argument("--spacing_maxvalue", "-smax", type = int, default= 30)
parser.add_argument("--spacing_minvalue", "-smin", type = int, default= 0)
parser.add_argument("--offset_left", "-ol", type = int, default=0)
parser.add_argument("--offset_right", "-or", type = int, default=0)
parser.add_argument("--len_motif", "-lm", type = int, default=10)
parser.add_argument("--output", "-o", type = str, default="interdistances")
parser.add_argument("--write_inter", "-wi",action='store_true', default= False)
parser.add_argument("--no_absolute_panel", "-nap",action='store_true', default= False)

args = parser.parse_args()

#print(args.negative_sets) #Debug print
main(	args.positive,
		args.negative_sets,
		args.threshold,
		args.spacing_maxvalue,
		args.spacing_minvalue,
		args.offset_left,
		args.offset_right,
		args.len_motif,
		args.output,
		args.write_inter,
		args.no_absolute_panel
	)
	
#python get_interdistances.py -pos /home/304.6-RDF/Jeremy/Arf_Paper_fig1_maker/results/2018_03_29_ARF2_fullset/ARF2_MP/Find_motifs/ARF2/scores/testing_pos_set.fasta.scores -neg /home/304.6-RDF/Jeremy/Arf_Paper_fig1_maker/results/2018_03_29_ARF2_fullset/ARF2_MP/Find_motifs/ARF2/scores/testing_neg_set.fasta.scores -ol 2 -or 2


