# -*- coding: utf-8 -*-
#!usr/bin/python

#################################################################################
#description     :This script will calculate the scores and density from a fasta file with a matrix file
#author          :Laura Gregoire
#date            :20160119
#version         :0.1
#usage           :python 4-scores.py -f <fasta file> -m <matrix file> 
#python_version  :2.7.3
# INPUT           :a fasta file, and the matrix in Scripts/files/matrix.txt
# OUTPUT          :2 files .fasta.scores and .fasta.dens
#################################################################################

import argparse
from argparse import RawTextHelpFormatter
import re
import numpy as np
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from os import popen
from math import *
from os.path import basename
from os.path import join

############
def analysisMatrix(matTxt): 
	""" to check the validity of the matrix
	# INPUT : matrix file content
	# OUTPUT : numbers according to the information of the matrix file"""
	first = matTxt[0].split()
	dep = 0
	if len(first)==5:
		name = first[3]
		if first[0] == "MATRIX":
			dep += 1000 
		if first[2] == "SYMMETRIC":
			dep +=10
		elif first[2] == "ASYMMETRIC":
			dep +=20
		if first[4] == "DEPENDENCY":
			dep += 1
		if first[1] == "SCORE":
			dep += 100 
		elif first[1] == "FREQUENCY":
			dep += 200 
		elif first[1] == "COUNT":
			dep += 300
		else:
			dep = 0
	second = matTxt[1].replace('\t', '')
	if second != "ACGT" :
		dep=0
	return dep,name

############
def error_matrix():
	""" error message if the matrix is not conform """
	print ""
	print "First line"
	print "~~~~~~~~~~~~~~"
	print "MATRIX\t[SCORE,FREQUENCY,COUNT]\t[SYMETRIC,ASYMMETRIC]\t{NAME}\t[DEPENDENCY/SIMPLE]"
	print "- MATRIX : needed"
	print "- [SCORE,FREQUENCY,COUNT] : to inform which kind of matrix this is"
	print "- [SYMETRIC,ASYMMETRIC] : to inform if the matrix is symetrical or not"
	print "- {NAME} : name of the matrix"
	print "- [DEPENDENCY/SIMPLE] : to inform if there are dependencies within the matrix"
	print ""
	print "Second line"
	print "~~~~~~~~~~~~~~"
	print "nucleotides separated by a tabulation (\\t):"
	print "A  C  G  T"
	print ""
	print "If dependencies:"
	print "~~~~~~~~~~~~~~"
	print "DEPENDENCY ([list of dependent nucleotides])"
	print "AA"
	print "AC"
	print "..."
	print "keeping the nucleotides order from line 2"

############
def creationtuplesreference(): 
	""" this function creates a tuple with all triplet in order 
	# INPUT : NO
	# OUTPUT : all the nucleotides arrangments for size = 1,2 et 3 """
	listebase = ['A','C','G','T']
	listedoublet = []
	listetriplet = []
	for base1 in listebase:
		for base2 in listebase:
			listedoublet.append(base1 + base2)
			for base3 in listebase:
				listetriplet.append(base1 + base2 + base3)
	return  tuple(listebase), tuple(listedoublet), tuple(listetriplet)
	
	
############
def freqToScore(matF,pseudoCount=1.0/600) :
	""" transform a frequency matrix into a score matrix
	Sc = ln(f+pC/fmax+pC)
	# INPUT: frequency matrix
	# OUTPUT: scoring matrix """
	matS = []
	for line in matF:
		lineFreq = []
		fmax = float(max(line)) + pseudoCount
		for freq in line:
			freq += pseudoCount
			tmp = np.log(float(freq)/fmax)
			lineFreq.append(tmp)
		matS.append(lineFreq)
	return matS
	
############
def countToScore(matC) :
	""" transform a counting matrix into a score matrix
	Sc = ln(f+pC/fmax+pC)
	# INPUT: counting matrix
	# OUTPUT: scoring matrix """
	matF = []
	pseudo = 0
	usePseudo = False
	for line in matC:
		lineCount = []
		counts = sum(line)
		pseudo = max(pseudo,counts)
		for count in line:
			tmp = count/counts
			lineCount.append(tmp)
			if count == 0:
				usePseudo = True
		matF.append(lineCount)
	if usePseudo:
		matS = freqToScore(matF,1/pseudo)
	else:
		matS = freqToScore(matF, 0)
	return matS

############
def modifMatrix(mat,FEAT,name) :	
	""" function to analyse and convert the matrix for the analysis of the data"""
	tuplebase, tupledoublet, tupletriplet = creationtuplesreference()
	matRes = []
	print "  Transcription Factor analysed %s"%(name)
	matRes.append(name)
	### remove empty lines
	tmp = []
	for line in mat:
		if line.replace('\t','') != 'ACGT':
			if line != '' :
				tmp.append(line)
	mat = tmp
	del(tmp)
	### delete useless lines
	mat.pop(0) # first line (with description)
	if "CORRELATION" in mat[-1]:# delete the first line (CORRELATION) if existes
		mat.pop(-1) 
	### counting the dependencies
	nDep = 0
	if int(FEAT[3]) == 1:
		for line in mat :
			if "DEPENDENCY" in line :
				nDep += 1
	print "  There is %s dependencies within the matrix"%(nDep)
	### taille du motif
	size = 0
	if nDep > 0 :
		while "DEPENDENCY" not in mat[size]:
			size += 1
	else:
		size = len(mat)
	print "  The binding motif is %s nucleotides long"%(size)
	matRes.append(size)
	matRes.append(1) # indication that the matrix is conform
	### keep the general matrix
	matGen = mat[:size]
	tmp = []
	for line in matGen:
		lineTmp = []
		for i in line.split():
			lineTmp.append(float(i))
		tmp.append(lineTmp)
	matGen = tmp
	del(tmp)
	### converting of the matrix if needed
	if int(FEAT[1]) == 3 :
		pwm = countToScore(matGen)
	elif int(FEAT[1]) == 2 :
		pwm = freqToScore(matGen)
	elif int(FEAT[1]) == 1 :
		pwm=matGen
	matRes.append(pwm)
	### keep the dependencies
	matDep = []
	if nDep > 0 :
		i=1
		deb = 0
		fin = 0
		while i <= nDep:
			if deb == 0 :
				deb = (size*i)-(2*(i))+2
			else :
				deb = fin+1
			fin = deb+16
			matDep.append(mat[deb:fin+1])
			i += 1
		del(i)
		matTmp = []
		for linesdep in matDep:
			dep = []
			loc = ""
			mats = []
			for line in linesdep:
				# keep the information about the location of the dependencies
				if "DEPENDENCY" in line:
					tmp = line.replace('DEPENDENCY (', '')
					loc = tmp.replace(')','')
					dep.append(loc.split())
				# analysis of the dependency matrices
				else:
					tmp = line.split()
					if re.search('[ACGT]',tmp[0]):
						tmp.pop(0)
					mats.append(tmp)
			dep.append(mats)
			# modification analysis of the dependency matrices
			if int(FEAT[1]) == 3 :
				pwm = countToScore(dep[1])
			elif int(FEAT[1]) == 2 :
				pwm = freqToScore(dep[1])
			elif int(FEAT[1]) == 1 :
				pwm=dep[1]
			tmp = {}
			if len(pwm) == 16:# three nucleotids dependency
				for i in range(0,len(pwm)):
					for j in [0,1,2,3]:
						tmp[tupletriplet[4*i+j]] = pwm[i][j]
			elif len(pwm) == 4: # two nucleotids dependency
				for i in range(0,len(pwm)):
					for j in [0,1,2,3]:
						tmp[tupledoublet[4*i+j]] = pwm[i][j]
			dep[1]=tmp
			matTmp.append(dep)
		# keep all the information in one variable
		matRes.extend(matTmp)
	return matRes
	
############
def SeqToDensities(fasta,mat,sym,output,matname,name="",chrom="",th=""):
	""" calculated the scores at all positions and the density of a sequence and retrieves those two informations within 2 separate files 
	# INPUT: fasta file, matrix, symmetry of the matrix and penalty
	# OUTPUT: scores and densities file"""
	sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        if not matname:
                fileScores = open(join(output,basename("%s.scores"%fasta)),"w")
        else:
                fileScores = open(join(output,name+"_"+basename("%s.scores"%fasta)),"w")
#	fileScores.write("name\tposition\tsens\tscore\n")
	# fileDensity = open(join(output,basename("%s.dens"%fasta)),"w")
	# fileDensity.write("name\tth\tdensity\n")
	print "  There is %s sequences to analyze"%(len(sequences))
	for s in sequences:
		seq=sequences[s].seq
		if len(chrom) != 0:
			if s.find(chrom) == -1:
				continue
		lSeq=len(seq)
		scores = SeqToScores(seq,mat)
		i=0
		for sc in scores:
			i=i+1
			if len(th) != 0:
				if sc <= float(th) or sc == "NA":
					continue
			fileScores.write("%s\t%s\t1\t%s\n"%(s,i,sc))
		if int(sym) == 2: # if the matrix is assymetric, look at both strands
#			print(mat)
#			print("\n")
			matRev = EvertedMatrix(mat)
			# print(matRev)
			scoresrev = SeqToScores(seq,matRev)
			scores.extend(scoresrev)
			i=0
			for sc in scoresrev:
				i=i+1
				if len(th) != 0:
					if sc <= float(th) or sc == "NA":
						continue
				fileScores.write("%s\t%s\t-1\t%s\n"%(s,i,sc))
		# keep the density for each region
		maxi = scoreMax(scores)
		scores = np.array(scores)
	# 	if maxi != "NA":
	# 		for i in range(long(maxi),0):
	# 			subset = scores[scores[:]>=i]
	# 			subset = subset[subset[:]<i+1]
	# 			n=len(subset)
	# 			fileDensity.write("%s\t%s\t%s\n"%(s,i,round(float(n)/len(seq)*1000,2)))
	# fileDensity.close()
	fileScores.close()

def returnMatrix(matrix) :
	matTmp = []
	for i in range(len(matrix)-1,-1,-1):
		line = matrix[i]
		lineT = [line[3],line[2],line[1],line[0]]
		matTmp.append(lineT)
	return matTmp


############
def EvertedMatrix(mat):	
	matRev = [mat[0],mat[1],mat[2],returnMatrix(mat[3])]
	if len(mat)>4:
		# modif position dependency
		for i in range(4,len(mat)):
			key = [int(j) for j in mat[i][0]]
			key = mat[1]-np.array(key)+1
			key = key.tolist()
			key = sorted(key)
			# modif val dependency
			tmp = mat[i][1]
			tmp2 = {}
			for keyTmpL in tmp.keys():
				keyTmp = list(keyTmpL)
				key2 = []
				for t in keyTmp:
					if t == 'A':
						key2.append("T")
					if t == 'C' :
						key2.append("G")
					if t == 'G' :
						key2.append("C")
					elif t == 'T' :
						key2.append("A")
				key2 = "".join(key2)[::-1]
				tmp2[key2] = tmp[keyTmpL]
			matRev.append([key,tmp2])
	return matRev

############

def SeqToScores(seq,mat) :
	""" calculates all the scores for one sequence
	# INPUT : sequence and the matrix
	# OUTPUT : list of scores"""
	tuplebase, tupledoublet, tupletriplet = creationtuplesreference()
	size, pwm= mat[1], mat[3]
	if len(seq)>=size :
		scores = []
		loc = []
		if len(mat)>4: #if dependencies keep the positions
			loc = []
			pwmD = mat[4:]
			for dep in pwmD:
				for i in dep[0]:
					j = long(i)-1
					loc.append(j)
		for c in range(len(seq)-size+1): # for all the positions within the sequence
			nucl = seq[c:c+size].upper() # subsequence of the motif size
			test = 0
			for nu in nucl :
				if nu not in ["A","C","G","T"]:
					test=1
			if test == 1:
				score = -70
			else :
				n=0
				score=0
				while n<size: 
					if (loc != []) and (n in loc): # score = dependencies score
						dep = pwmD[(loc.index(n)/len(pwmD))]
						lDep = dep[0]
						mDep = dep[1]
						subseq = nucl[(long(lDep[0])-1):long(lDep[-1])]
						score = score + float(mDep[str(subseq)])
						n += len(lDep)
					else:# score = dependencies score + score
						score = score + float(pwm[n][tuplebase.index(nucl[n])])
						n += 1
			scores.append(score)
	else: # in case seq is to small just add a very negative score so script doesn't bug.
		scores=[-66]
	return scores
		
############
def scoreMax(scores) :
	""" retrieves the manimum score within a list of scores"""  
	tmp = []
	for el in scores:
		if el != "NA" :
			tmp.append(el)
	if tmp == [] :
		return -70
	else:
		return floor(min(tmp))
##############################################################################
def main(fastaFile,matFile,output,matname,chrom,th):
	if fastaFile != None or matFile != None :
		# analysis of the matrix
		print "Analyse de la matrice"
		mat = open(matFile, "r").read()
		mat = mat.replace("\r","").split("\n")
		valMat,name = analysisMatrix(mat)
		if valMat == 0 :
			print "  La matrice n'est pas conforme et doit ressembler à ça :"
			error_matrix()
		else:
			print "  La matrice est conforme"
			FEAT = list(str(valMat))
			#FEAT[0] : MATRIX
			#FEAT[1] : SCORE/FREQUENCY/COUNT
			#FEAT[2] : SYMMETRIC/ASYMMETRIC
			#FEAT[3] : DEPENDENCY/SIMPLE
			matRes = modifMatrix(mat,FEAT,name)
		    # description of the 'matRes' resulting matrix:
		    # [name,size,conformity,[[.,.,.,.],...]general matrix,
		    # [[[loc],[[.,.,.,.],...],...]dependency matrix]
		    # len = 4 if no dependencies
		    # len > 4 if dependencies
			if FEAT[1] != "1":
				matExp = open(join(output,"PWM"+matRes[0]+".txt"),"w")
				matNew=matRes[3]
				desc = mat[0].split()
				matExp.write(desc[0]+" SCORE "+desc[2]+" "+desc[3]+" "+desc[4]+"\n")
				matExp.write("A\tC\tG\tT\n")
				for line in matNew:
					exp = []
					for el in line:
						exp.append(str(el))
					matExp.write("\t".join(exp)+"\n")
				if len(matRes)>4:
					for d in range(4,len(matRes)):
						matExp.write("\n")
						DepDict = {}
						matDep = matRes[d]
						pos = matDep[0]
						matExp.write("DEPENDENCY ("+" ".join(pos)+")\n")
						matExp.write("A\tC\tG\tT\n")
						for nucl in matDep[1]:
							First = nucl[:2]
							Sec = tuplebase.index(nucl[-1])
							if First not in DepDict :
								DepDict[First] = ["","","",""]
							tmp = DepDict[First]
							tmp[Sec] = str(matDep[1][nucl])
							DepDict[First] = tmp

						for d in tupledoublet:
							matExp.write(d+" "+" ".join(DepDict[d])+"\n")
				
				matExp.close()
			
	
		    ### calculation of scores and densities
			print "  Calculation of scores and densities"
                        SeqToDensities(fastaFile,matRes,FEAT[2],output,matname,name,chrom,th)
                                
	else:
		parser.print_help()

##############################################################################

if __name__ == "__main__": 
	print "\n      TRANSCRIPTION FACTOR ANALYSIS\n"
	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-f","--fasta", type=str, help='give the fasta file with the sequences')
	parser.add_argument("-m","--matrixfile", type=str, help='give the matrix file')
	parser.add_argument("-o","--output", type=str, help='give output file path and name')
        parser.add_argument( "-matname","--matrix_name",action='store_true', default= False)
	parser.add_argument("-chr","--chromosome", type=str, help='give chr to conserve',default="")
        parser.add_argument("-th", "--threshold",type=str, help='give threshold to remove bad ones',default="")

	args = parser.parse_args()
	
	fastaFile = args.fasta
        matname = args.matrix_name
	matFile = args.matrixfile
	output = args.output
	chrom = args.chromosome
	th = args.threshold
	main(fastaFile,matFile,output,matname,chrom,th)

