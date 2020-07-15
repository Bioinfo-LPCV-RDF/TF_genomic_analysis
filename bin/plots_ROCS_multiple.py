# -*- coding: utf-8 -*-
#!usr/bin/python

################################################                General intels                ################################################                
#description     : 
#author          :Jérémy LUCAS inspired by a work from Arnaud STIGLIANI
#date            : 17-01-2018
#version         :0.1
#usage           : 
#python_version  :2.7.9
################################################                import packages                ################################################                

import argparse
import os
import matplotlib.pyplot as plt
import sklearn.metrics
import numpy as np
################################################                Parsing arguments                ################################################

parser = argparse.ArgumentParser() 

parser.add_argument("--scores","-s" , nargs="+")
parser.add_argument("--names","-n" , nargs="+")
parser.add_argument("--colors","-c" , nargs="+")
parser.add_argument("--output","-o")
parser.add_argument("--outputfile","-of")
args = parser.parse_args()
name_scores_files = args.scores
names = args.names
dir_OUT = args.output
fileout = args.outputfile
colors = args.colors
################################################                Functions                ################################################'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black']
def prepare_curve(groups,scores):
    fpr,tpr,_=sklearn.metrics.roc_curve(groups,scores, pos_label=1)
    auc=sklearn.metrics.auc(fpr,tpr)
    auc=round(auc,3)
    return fpr,tpr,auc

def plot_curve(scores,name_scores_files,names,dir_OUT,fileout,colors):
    fig=plt.figure(figsize=(6, 6), dpi=60)
    ax=fig.add_subplot(111)
    ax.set_xlabel("fraction of unbound sequences", fontsize=14)
    ax.set_ylabel("fraction of bound sequences", fontsize=14)
    i=0
    #colors=['#40A5C7', '#307C95','#F9626E', '#BB4A52']#["blue","green","red","black","purple","orange"]
	# SEP3 '#F0875A', '#B46544'
	# SEP3AG '#40A5C7', '#307C95'
	#SEP3delAG '#F9626E', '#BB4A52'
    for elt in name_scores_files:
    	ax.plot(scores[elt]["fpr"],scores[elt]["tpr"], colors[i], linewidth=3)
    	ax.annotate("AUC "+names[i]+" = "+str(scores[elt]["auc"]),xy=(0.2,0.2),xytext=(0.15,0.05+(0.05*i)), color=colors[i], fontsize=18)
    	i+=1
#    ax.plot(fpr_2,tpr_2, 'blue')
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle='dashed', color='black')
    ax.axis([0.0, 1.0,0.0, 1.0])

   
    
#    ax.annotate("AUC "+name_2+" = "+str(auc_2),xy=(0.2,0.2),xytext=(0.4,0.2), color='blue')
    
    plt.tight_layout()
#    plt.show()
    plt.savefig(dir_OUT+fileout)
################################################                MAIN                ################################################

scores={}
fpr=[]
tpr=[]
auc=[]
for elt in name_scores_files:
	scores[elt]={}
	scores[elt]["groups"]=[]
	scores[elt]["scores"]=[]
	scores[elt]["fpr"]=[]
	scores[elt]["tpr"]=[]
	scores[elt]["auc"]=0
	
	
	with open(elt,"r") as IN:
		for line in IN:
			columns=line.split("\t")
			if columns[0] != "NA" and columns[1] != "NA":
				scores[elt]["groups"].append(1)
				scores[elt]["scores"].append(float(columns[0]))
				scores[elt]["groups"].append(2)
				scores[elt]["scores"].append(float(columns[1]))
	(scores[elt]["fpr"],scores[elt]["tpr"],scores[elt]["auc"])=prepare_curve(scores[elt]["groups"],scores[elt]["scores"])
	
plot_curve(scores,name_scores_files,names,dir_OUT,fileout,colors)


#scores_1=[]
#groups_1=[]
#fpr_1=""
#tpr_1=""
#auc_1=""

#IN.close()

#scores_2=[]
#groups_2=[]
#fpr_2=""
#tpr_2=""
#auc_2=""
#with open(name_scores_file_2,"r") as IN:
#    for line in IN:
#        columns=line.split("\t")
#        if columns[0] != "NA" and columns[1] != "NA":
#            groups_2.append(1)
#            scores_2.append(float(columns[0]))
#            groups_2.append(2)
#            scores_2.append(float(columns[1]))
#IN.close()

#(fpr_1,tpr_1,auc_1)=prepare_curve(groups_1,scores_1)
#(fpr_2,tpr_2,auc_2)=prepare_curve(groups_2,scores_2)
#plot_curve(fpr_1,tpr_1,auc_1,fpr_2,tpr_2,auc_2,name_1,name_2)

