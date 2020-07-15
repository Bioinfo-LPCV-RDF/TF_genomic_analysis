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
import numpy as np
import string
from Bio import SeqIO # version 1.70

################################################                Functions                ################################################

################################################                MAIN                ################################################
def main(gffin,gffout,fastafile):
	table_types={}
	gff=np.genfromtxt(gffin,dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')]) # read gff
	types=np.unique(np.sort(gff['type'])) # to get each type

	if 'chromosome' not in types:
		lines=np.empty(0,dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')])
		fasta = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
		for elt in fasta.keys():
			length=len(fasta[elt].seq)
			line= np.array([(str(elt),"LPCV","chromosome",1,length,".",".",".","ID="+str(elt)+";Name="+str(elt))],dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')])
			lines=np.concatenate((lines,line))
		gff = np.concatenate((gff,lines))
	gff=np.sort(gff,order=['chrom','start','stop'])
	np.savetxt(gffout,gff,delimiter="\t",fmt="%s")
	
################################################                Parsing arguments                ################################################

parser = argparse.ArgumentParser() 

parser.add_argument("--gffin","-gi",type=str)
parser.add_argument("--gffout","-go",type=str)
parser.add_argument("--fasta","-f",type=str)

args = parser.parse_args()
main(args.gffin,args.gffout,args.fasta)
