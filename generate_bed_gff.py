#-*- coding: utf-8 -*-
import numpy as np
import numpy.core.defchararray as np_f
import argparse
from argparse import RawTextHelpFormatter
import os

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-g","--gff", type=str, help='give the entry gff file')
parser.add_argument("-p","--promoter", type=int, help='promoter size')
parser.add_argument("-o","--output", type=str, help='give output file directory and name')
args = parser.parse_args()

gff = args.gff
promoter_size = args.promoter
output = args.output


table_types={} #to clssify regions by type
gff=np.genfromtxt(gff,dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')]) # read gff
############# to remove alternative splicings ##########
# if you want to keep alternatively spliced rna, you can
# comment the 4 following lines
Ind=np.char.find(gff['ID'],'.')
Ind2=np.char.find(gff['ID'],'.1')
Ind3=np.where(np.logical_and(Ind != -1,Ind2 == -1))
gff=np.delete(gff,Ind3)
#######################################################
types=np.unique(np.sort(gff['type'])) # to get each type
chromosomes=gff[np.where(gff['type']=="chromosome")] # get chromosome informations
all_gff=np.empty(0,dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')])#initialize final gff

  
for elt in chromosomes:
    chrom=elt[0]
    begin=elt[3]
    end=elt[4]
    tab=gff[np.where(gff['chrom']==chrom)]
    for elt2 in types:
        table_types[elt2]=tab[np.where(tab['type']==elt2)] # retrieve regions from gff and store them in table_type
    ############### Promoters ######################
    promoter_plus=table_types['gene'][np.where(table_types['gene']['strand']=='+')] # define + strand promoters
    promoter_minus=table_types['gene'][np.where(table_types['gene']['strand']=='-')] # define - strand promoters
    promoter_plus['stop']=promoter_plus['start'] # if + strand, promoter stop = begin of gene
    promoter_plus['start']=promoter_plus['start'] - promoter_size # promoter start = begin of gene - promoter size
    promoter_plus['start'][np.where(promoter_plus['start'] < 0)] = 0 # be careful about chromosome edges
    promoter_minus['start']=promoter_minus['stop'] # same but for other strand
    promoter_minus['stop']=promoter_minus['stop'] + promoter_size
    promoter_plus['stop'][np.where(promoter_plus['stop'] > table_types['chromosome']['stop'][0])] = table_types['chromosome']['stop'][0]
    promoters=np.concatenate((promoter_minus,promoter_plus),0) # gather + and - promoter 
    promoters['type']='promoter' # create promoter type
    table_types['promoter']=promoters # add promoters type to table_type
    ###############################################
    ############### Introns #######################
    if 'exon' not in types and 'CDS' in types:
        table_types['exon']=np.copy(table_types['CDS'])
    tab1=np.delete(table_types['exon'],[len(table_types['exon'])-1]) # create exon table w/o last element
    tab2=np.delete(table_types['exon'],[0]) # create exon table w/o first element
    a=np.array(np.char.split(tab1['ID'],sep="Parent=").tolist())
    b=np.hsplit(a,2)[1]
    tab1['ID']=b.flatten()
    a=np.array(np.char.split(tab2['ID'],sep="Parent=").tolist())
    b=np.hsplit(a,2)[1]
    tab2['ID']=b.flatten()
    Ind_pos=np.where(np.logical_and(tab1['ID']==tab2['ID'], np.logical_and(tab1['stop'] < tab2['start'],tab1['strand']=='+') )) # Compare start and stop of succesive exons of same gene  
    introns_pos=tab1[Ind_pos] # retrieve exon when there is a gap between two exons 
    introns_pos['start']=introns_pos['stop'] # adjust boundaries to create corresponding intron
    introns_pos['stop']=tab2[Ind_pos]['start'] # adjust boundaries to create corresponding intron
    introns_pos['type']='intron' # create intron type
    
    Ind_neg=np.where(np.logical_and(tab1['ID']==tab2['ID'], np.logical_and(tab2['stop'] < tab1['start'],tab1['strand']=='-') )) # Compare start and stop of succesive exons of same gene  
    introns_neg=tab1[Ind_neg] # retrieve exon when there is a gap between two exons 
    introns_neg['stop']=introns_neg['start'] # adjust boundaries to create corresponding intron
    introns_neg['start']=tab2[Ind_neg]['stop'] # adjust boundaries to create corresponding intron
    introns_neg['type']='intron' # create intron type
#    print(introns_neg)
#    print(introns_pos)
    introns=np.hstack((introns_pos,introns_neg))
    table_types['intron']=introns # add intron type to table_type
    ###############################################
    del(table_types['gene']) # del gene type
    if 'protein' in types:
        del(table_types['protein']) # del gene type
    if 'mRNA' in types:
        del(table_types['mRNA']) # del gene type
    if 'CDS' in types:
        del(table_types['CDS'])# del CDS type
    del(table_types['chromosome']) # del chromosome type
    table_types['5UTR']=np.copy(table_types['five_prime_UTR']) # change UTR keys 
    table_types['3UTR']=np.copy(table_types['three_prime_UTR'])
    del(table_types['five_prime_UTR'])
    del(table_types['three_prime_UTR'])
    table_types['5UTR']['type']='5UTR'
    table_types['3UTR']['type']='3UTR'
    all_gff_chrom=np.empty(0,dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')]) # create final gff for corresponding chromosome
    for origine in table_types:
        if (len(table_types[origine]) != 0): 
            all_gff_chrom=np.concatenate((all_gff_chrom,table_types[origine])) # add each table_type entry to all_gff_chrom if there is elt in entry
    ############### Intergenic #######################
    all_gff_chrom=np.sort(all_gff_chrom,order=['start','stop']) # sort file before looking for gaps (intergenics)
    tab2=np.delete(all_gff_chrom,[0]) # same method as for intron research
    tab1=np.delete(all_gff_chrom,[len(all_gff_chrom)-1])
    where_intergenics=np.where(tab2['start']>tab1['stop'])
    all_intergenics=np.zeros(len(where_intergenics[0]),dtype=[('chrom','S20'),('genome','S300'),('type','S300'),('start',int),('stop',int),('col1','S5'),('strand','S5'),('col2','S5'),('ID','S300')])
    all_intergenics['chrom']=chrom
    all_intergenics['type']='intergenic'
    all_intergenics['start']=tab1[where_intergenics]['stop']
    all_intergenics['stop']=tab2[where_intergenics]['start']
    ###################################################
    all_gff_chrom=np.concatenate((all_gff_chrom,all_intergenics)) # add intergenics to all_gff_chrom
    all_gff=np.concatenate((all_gff,all_gff_chrom)) # add all_gff_chrom to final gff

#all_gff['chrom']=np_f.replace(all_gff['chrom'],'Chr','chr') # change Chr to chr
bed=np.zeros(len(all_gff),dtype=[('chrom','S20'),('start',int),('stop',int),('type','S300')]) # create bed file with needed informations
bed['chrom']=all_gff['chrom']
bed['type']=all_gff['type']
bed['start']=all_gff['start']
bed['stop']=all_gff['stop']
np.savetxt(output,bed,delimiter="\t",fmt="%s") # save bed file in output directory



















