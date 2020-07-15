import re
import numpy as np
from Bio import SeqIO

####################### This is useful to calculate interdependances ###########

codigo = { 'ACC': 5, 'ATG': 14, 'AAG': 2, 'AAA': 0, 'ATC': 13, 'AAC': 1, 'ATA': 12,
           'AGG': 10, 'CCT': 23, 'ACT': 7, 'AGC': 9, 'ACA': 4, 'AGA': 8, 'CAT': 19, 
           'AAT': 3, 'ATT': 15, 'CTG': 30, 'CTA': 28, 'CTC': 29, 'CAC': 17, 'ACG': 6,
           'CAA': 16, 'AGT': 11, 'CCA': 20, 'CCG': 22, 'CCC': 21, 'TAT': 51, 'GGT': 43, 
           'TGT': 59, 'CGA': 24, 'CAG': 18, 'CGC': 25, 'GAT': 35, 'CGG': 26, 'CTT': 31, 
           'TGC': 57, 'GGG': 42, 'TAG': 50, 'GGA': 40, 'TAA': 48, 'GGC': 41, 'TAC': 49,
           'GAG': 34, 'TCG': 54, 'TTA': 60, 'GAC': 33, 'CGT': 27, 'TTT': 63,
           'TCA': 52, 'GCA': 36, 'GTA': 44, 'GCC': 37, 'GTC': 45, 'GCG': 38,
           'GTG': 46, 'TTC': 61, 'GTT': 47, 'GCT': 39, 'TGA': 56, 'TTG': 62,
           'TCC': 53, 'TGG': 58, 'GAA': 32, 'TCT': 55}
           
codigoi = { "A" : "T", "C" : "G", "G" : "C", "T" : "A"}

###############################################################################

def divide(a, b):
	if b == 0:
		return np.nan
	else: 
		return a/b

def get_score_matrix(Mdata,matrixType,pseudoCount) :
	if matrixType == "freq" :
		## These lines allows to transform the frequency values into scores values
		matScore = []
		lenMotif = 0
		a = 0
		for i in range(0,len(Mdata)/4):
			lenMotif = lenMotif + 1
			fmax = float(max(float(Mdata[a]),float(Mdata[a+1]),float(Mdata[a+2]),float(Mdata[a+3]))) #+ pseudoCount
			for j in range (0,4):
				matScore.append(np.log(float(float(Mdata[a+j]) ) /fmax)) # can add pseudoCount to numerator
			a = a + 4
	else :
		matScore = map(float,Mdata)
		lenMotif = len(matScore)/4
	print(matScore)
	return (matScore, lenMotif)

def get_dependency_matrix(dependencyFile,num) : 
	G = open(dependencyFile,"r")
	dependency_file_content = G.read().replace("\r","\n") + "\n"
	G.close()
	
	num2 = re.compile(r"(?<![\d.])[0-9]+(?![\d.])")
	position_dependency_matrix = num2.findall(dependency_file_content)
	position_dependency_matrix = map(int, position_dependency_matrix)
	
	dependency_matrix = num.findall(dependency_file_content)
	dependency_matrix = map(float, dependency_matrix)
	
	dependency_matrix_associated_with_position = []
	index1 = 0
	index2 = 3
	index3 = 0
	index4 = 64
	for i in range(0, 3):
		dependency_matrix_associated_with_position.append(position_dependency_matrix[index1:index2])
		dependency_matrix_associated_with_position.append(dependency_matrix[index3:index4])
		index1 = index1 + 3
		index2 = index2 + 3
		index3 = index3 + 64
		index4 = index4 + 64
		
	return(dependency_matrix_associated_with_position)

def seq_c(site):
        site_i = site[-1::-1]
        site_c = ""
        for x in site_i:
		y = codigoi[x]
                site_c = site_c + y
        return site_c   

def add_scores_associated_with_interdependent_positions(dependency_matrix,scoreStrandPos,scoreStrandNeg,strandPos):
	cStrand = ""
	for lettre in seq_c(strandPos):
		cStrand = lettre + cStrand
	cStrand = cStrand[::-1]
	
	site1 = ""
	Csite1 = ""
	site2 = ""
	Csite2 = ""
	site3 = ""
	Csite3 = ""
	for i in dependency_matrix[0]:
		site1 = site1 + strandPos[i-1]
		Csite1 = Csite1 + cStrand[i-1]
	for i in dependency_matrix[2]:
		site2 = site2 + strandPos[i-1]
		Csite2 = Csite2 + cStrand[i-1]
	for i in dependency_matrix[4]:
		site3 = site3 + strandPos[i-1]
		Csite3 = Csite3 + cStrand[i-1]
	#print("dependency_matrix[1][0] : ",dependency_matrix[1][0])
	scoreStrandPos = scoreStrandPos + dependency_matrix[1][codigo[site1]] + dependency_matrix[3][codigo[site2]] + dependency_matrix[5][codigo[site3]] 
	scoreStrandNeg = scoreStrandNeg + dependency_matrix[1][codigo[Csite1]] + dependency_matrix[3][codigo[Csite2]] + dependency_matrix[5][codigo[Csite3]]
	return(scoreStrandPos, scoreStrandNeg)

def interdistance_calcul(InterDR,InterER,InterIR,sum_threshold,good_score_positions,offset_left,offset_right,lenMotif,Interdistance_maxValue,write_inter,f,threshold) :
        list_d=[]
	for first in range(0,len(good_score_positions)-1) :
		firstSubSeq = good_score_positions[first]
		for second in range(first+1,len(good_score_positions)) :
    			        secondSubSeq = good_score_positions[second]

				'''
				Here we do different substractions according to we get a Direct Repeat (DR), an Inverted Repeat (IR) or an Everted Repeat (ER).
				Because In the litterature, The interdistance calculation was made between subsequences from our matrix
				and not between the whole sequences from our matrix.
				So according to your specific interdistance calculations you can change these lines.
				'''

				if sum_threshold == True :
					if int(firstSubSeq[2]) + int(secondSubSeq[2]) > sum_thresold :
						if firstSubSeq[1] == ">" and secondSubSeq[1] == ">" :
							d = ( int(secondSubSeq[0]) +offset_left ) -( int(firstSubSeq[0]) + lenMotif -offset_right )
							if Interdistance_maxValue >= d >= 0 :
								InterDR[d] += 1
						if firstSubSeq[1] == "<" and secondSubSeq[1] == "<" :
							d = ( int(secondSubSeq[0]) +offset_right ) -( int(firstSubSeq[0]) + lenMotif -offset_left)
							if Interdistance_maxValue >= d >= 0 :
								InterDR[d] += 1
						if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
							d = ( int(secondSubSeq[0]) +offset_right ) -( int(firstSubSeq[0]) + lenMotif -offset_right )
							if Interdistance_maxValue >= d >= 0 :
								InterER[d] += 1
						if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
							d = ( int(secondSubSeq[0]) +offset_left ) -( int(firstSubSeq[0]) + lenMotif -offset_left )
							if Interdistance_maxValue >= d >= 0 :
								InterIR[d] += 1
				else :
					if firstSubSeq[1] == ">" and secondSubSeq[1] == ">" :
						d = ( int(secondSubSeq[0]) +offset_left ) -( int(firstSubSeq[0]) + lenMotif -offset_right )
						if Interdistance_maxValue >= d >= 0 :
							InterDR[d] += 1
                                                        if write_inter is True :
                                                                list_d.append("DR"+str(d)+"/"+str(firstSubSeq[2])+":"+str(secondSubSeq[2]))#+"/"+str(firstSubSeq[0])+":"+str(secondSubSeq[0]))
					if firstSubSeq[1] == "<" and secondSubSeq[1] == "<" :
						d = ( int(secondSubSeq[0]) +offset_right ) -( int(firstSubSeq[0]) + lenMotif -offset_left)
						if Interdistance_maxValue >= d >= 0 :
							InterDR[d] += 1
                                                        if write_inter is True :
                                                                list_d.append("ER"+str(d)+"/"+str(firstSubSeq[2])+":"+str(secondSubSeq[2]))#+"/"+str(firstSubSeq[0])+":"+str(secondSubSeq[0]))
					if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
						d = ( int(secondSubSeq[0]) +offset_right ) -( int(firstSubSeq[0]) + lenMotif -offset_right )
						if Interdistance_maxValue >= d >= 0 :
							InterER[d] += 1
                                                        if write_inter is True :
                                                                list_d.append("IR"+str(d)+"/"+str(firstSubSeq[2])+":"+str(secondSubSeq[2]))#+"/"+str(firstSubSeq[0])+":"+str(secondSubSeq[0]))
					if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" : # 
						d = ( int(secondSubSeq[0]) +offset_left ) -( int(firstSubSeq[0]) + lenMotif -offset_left )
						if Interdistance_maxValue >= d >= 0 :
							InterIR[d] += 1
                                                        if write_inter is True :
                                                                list_d.append("DR"+str(d)+"/"+str(firstSubSeq[2])+":"+str(secondSubSeq[2]))#+"/"+str(firstSubSeq[0])+":"+str(secondSubSeq[0]))
        if write_inter is True :
                try:
                        if len(list_d)==0:
                                list_d=["distance"]
                        f.write(firstSubSeq[3].replace(":","\t").replace("-","\t")+"\t"+','.join(list_d)+"\n")
                except:
                        f.write("chr\t1\t1\tdistances\n")
	return(InterDR,InterER,InterIR)
	
def get_interdist(matF, matRev, FastaFile, threshold, offset_left,offset_right, Interdistance_maxValue, sum_threshold, lenMotif, dependencyFile, sequence_number,tffm_first_order,write_inter,output):
        if tffm_first_order != "":
                (ER,IR,DR)=get_interdist_tffm(FastaFile, threshold, offset_left,offset_right, Interdistance_maxValue, lenMotif, sequence_number,tffm_first_order,write_inter,output)
                return(ER,IR,DR)
	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))
        from Bio.SeqRecord import SeqRecord
	print "There are %s sequence(s) to analyze"%(sequence_number)
	# We will store in these lists all the occurences of each kind of interdistances between motifs found in all sequences.
	DR = [] 
	ER = [] 
	IR = []
        if write_inter is True:
                f= open(output+".bed","w")
                f.close()
                f= open(output+".bed","w+")                                
        else :
                f=""
	for a in threshold :
		DR.append([0] * (Interdistance_maxValue + 1) )
		ER.append([0] * (Interdistance_maxValue + 1) )
		IR.append([0] * (Interdistance_maxValue + 1) )
	
	score_occurence = 0
	nb = 0
	# We look at all the fasta sequences:
	for s in sequences:
		# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
		#if type(threshold) is list:
		good_score_positions = []
		ctrl_scores=[]
		if sum_threshold == False :
			for a in threshold :
				good_score_positions.append( [] )
		bestScore = 0
		positionOfTheBestScore = 0
		# This line allows to retrieve the DNA sequence
		seq = sequences[s].seq
		id=sequences[s].id
		# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
		for c in range(len(seq) - (lenMotif -1)):
			strandPos = seq[c:c+lenMotif].upper()
			test = 0
			for nu in strandPos :
				if nu not in ["A","C","G","T"]:
					test = 1
			if test == 1:
				score = "NA"
			else :
                                n = 0
                                #These lines allows to calculate a score for one sub-sequence
                                scoreStrandPos = 0
                                scoreStrandNeg = 0 
                                while n<lenMotif:
                                        if strandPos[n] == 'A':
                                                scoreStrandPos = scoreStrandPos + matF[n*4]
                                                scoreStrandNeg = scoreStrandNeg + matRev[n*4]
                                        elif strandPos[n] == 'C':
                                                scoreStrandPos = scoreStrandPos + matF[n*4+1]
                                                scoreStrandNeg = scoreStrandNeg + matRev[n*4+1]
                                        elif strandPos[n] == 'G':
                                                scoreStrandPos = scoreStrandPos + matF[n*4+2]
                                                scoreStrandNeg = scoreStrandNeg + matRev[n*4+2]
                                        elif strandPos[n] == 'T':
                                                scoreStrandPos = scoreStrandPos + matF[n*4+3]
                                                scoreStrandNeg = scoreStrandNeg + matRev[n*4+3]			
                                        n += 1
                                num = re.compile(r"([+-]?\d+[.,]\d+)")
                                if dependencyFile != "" : 			
                                        scoreStrandPos, scoreStrandNeg = add_scores_associated_with_interdependent_positions(get_dependency_matrix(dependencyFile,num),scoreStrandPos,scoreStrandNeg,strandPos)	
				#These lines allows to retrieve the position and the strand where there is a predicted binding site. 
				#You can change the threshold.
				if sum_threshold == True :
					if scoreStrandPos > min(threshold):
						good_score_positions.append([c+1,">",scoreStrandPos])
					if scoreStrandNeg > min(threshold):
						good_score_positions.append([c+1,"<",scoreStrandNeg])
				else :
					for a , b in zip(good_score_positions, threshold) :
						if scoreStrandPos > b:
							score_occurence = score_occurence + 1
							a.append([c+1,">",scoreStrandPos])
							if b == min(threshold):
								ctrl_scores.append(str(sequences[s].id)+"\t"+str(c+1)+"\t"+"1"+"\t"+str(scoreStrandPos))
						if scoreStrandNeg > b:
							score_occurence = score_occurence + 1
							a.append([c+1,"<",scoreStrandNeg])
							if b == min(threshold):
								ctrl_scores.append(str(sequences[s].id)+"\t"+str(c+1)+"\t"+"-1"+"\t"+str(scoreStrandNeg))
		with open(FastaFile+"scores_interdistv1.tsv","a") as OUT:
			for elt in ctrl_scores:
				OUT.write(str(elt)+"\n")
		# Once we have stored all the positions, we calculate all the interdistances:
		if sum_threshold == True :
			for interDIR, interEVER, interINVER,sum_threshold in zip(DR,ER,IR,threshold) :
				interdistance_calcul(interDIR,interEVER,interINVER,sum_thresold,good_score_positions,offset_left,offset_right,lenMotif,Interdistance_maxValue,False,"")
		else :
			for goodScores, interDIR, interEVER, interINVER,a in zip(good_score_positions,DR,ER,IR,threshold) :
				InterDR,InterER,InterIR = interdistance_calcul(interDIR,interEVER,interINVER,threshold,goodScores,offset_left,offset_right,lenMotif,Interdistance_maxValue,write_inter,f,a)
		if sequence_number :
			nb = nb + 1
		if nb == sequence_number : 
			break
        if write_inter is True:
                f.close()
	return(DR,ER,IR)

#def get_interdist_tffm(FastaFile, threshold, offset_left,offset_right, Interdistance_maxValue, lenMotif, sequence_number,tffm_first_order,write_inter,output):
#	DR = [] 
#	ER = [] 
#	IR = []
#	all_score_positions=[]
#	hit_list=[]
#	for a in threshold :
#		DR.append([0] * (Interdistance_maxValue + 1) )
#		ER.append([0] * (Interdistance_maxValue + 1) )
#		IR.append([0] * (Interdistance_maxValue + 1) )
#	saved=""
#	for hit_init in tffm_first_order.scan_sequences(FastaFile, only_best=False):
#		if hit_init:
#			hit=str(hit_init).strip("\n").split("\t")
#			if hit[0]!=saved:
#				saved=hit[0]
#				all_score_positions.append([])
#			if hit[3]=="+":
#				all_score_positions[-1].append([int(hit[1]),">",float(hit[7]),str(hit[0])])
#			elif hit[3]=="-":
#				all_score_positions[-1].append([int(hit[1]),"<",float(hit[7]),str(hit[0])])
#		
#	if write_inter is True:
#		f= open(output+".bed","w")
#		f.close()
#		f= open(output+".bed","w+")
#	else:
#		f=""
#	for seq  in all_score_positions:
#		good_score_positions=[]
#		for a in threshold :
#			good_score_positions.append([])
#			for score in seq :
#				if score[2] > a:
#					good_score_positions[-1].append(list(score))
#	for goodScores, interDIR, interEVER, interINVER, a in zip(good_score_positions,DR,ER,IR,threshold) :
#		InterDR,InterER,InterIR = interdistance_calcul(interDIR,interEVER,interINVER,threshold,goodScores,offset_left,offset_right,lenMotif,Interdistance_maxValue,write_inter,f,a)
#	if write_inter is True:
#		f.close()
#	return(DR,ER,IR)

def get_interdist_tffm(FastaFile, threshold, offset_left,offset_right, Interdistance_maxValue, lenMotif, sequence_number,tffm_first_order,write_inter,output):
	DR = [] 
	ER = [] 
	IR = []
	all_score_positions=[]
        hit_list=[]
	for a in threshold :
		DR.append([0] * (Interdistance_maxValue + 1) )
		ER.append([0] * (Interdistance_maxValue + 1) )
		IR.append([0] * (Interdistance_maxValue + 1) )
	for hit in tffm_first_order.scan_sequences(FastaFile, only_best=False):
                if hit:
                        hit_list.append(str(hit).strip("\n").split("\t"))
        saved=""
        for hit in hit_list:
                if str(hit[0])!=saved:
                        all_score_positions.append([])
                        saved=str(hit[0])
                if hit[3]=="+":
                        all_score_positions[-1].append([int(hit[1]),">",float(hit[7]),str(hit[0])])
                elif hit[3]=="-":
                        all_score_positions[-1].append([int(hit[1]),"<",float(hit[7]),str(hit[0])])
        if write_inter is True:
                f= open(output+".bed","w")
                f.close()
                f= open(output+".bed","w+")
        else:
                f=""                
        for seq  in all_score_positions:
                good_score_positions=[]
                for a in threshold :
                        good_score_positions.append([])
                        for score in seq :
                                if score[2] > a:
                                        good_score_positions[-1].append(list(score))
                for goodScores, interDIR, interEVER, interINVER, a in zip(good_score_positions,DR,ER,IR,threshold) :
			InterDR,InterER,InterIR = interdistance_calcul(interDIR,interEVER,interINVER,threshold,goodScores,offset_left,offset_right,lenMotif,Interdistance_maxValue,write_inter,f,a)
        if write_inter is True:
                f.close()
	return(DR,ER,IR)
