#coding: utf-8

from math import exp
import os
import argparse


#constants (according to Moyroud and al., 2011) (valid for SYM, full site matrix LFY.pfm)

#a=2.5663
#b=0.3598
#TF = exp(b/a)
#threshold = -20

#output = "/home/304.6-RDF/Pauline/POcc" 
#matrix = "/home/304.6-RDF/data/LFY.pfm"
#sequences = "/home/304.6-RDF/FParcy_prom_results/apetala3/promoters_all_clade/promoters_seq_by_length/prom1000.fa"

#voir pour passer le nom du gène en argument (facultatif, default = "")



def Ka(score,a,b) :
	return exp(-(b-float(score))/a)


	
def main(scores,output,a,b) : 
#	get_score_file="python /home/304.6-RDF/scripts/scores.py -f "+ sequences + " -m " + matrix + " -o " + output
#	os.system(get_score_file) # on a le fichier de scores
	
#	pocc_scores = open(output,"w") #create file if it doesn't already exist
#	pocc_scores.close() #clear file
	
	TF = exp(b/a)
	
	with open(scores, "r") as score_file, open(output,"w") as pocc_scores :
		s = 0
#		firstline = score_file.readline()
		current_gene = "" # firstline.split('\t')[0]
		
		for line in score_file :
			line = line.split('\t')
			if current_gene != "" and line[0] == current_gene : 
#				if float(line[3]) >= threshold :
				s+= Ka(line[3],a,b)*TF/(1+Ka(line[3],a,b)*TF) 
			else :
				pocc_scores.write(str(s).encode('utf-8')+"\n") #on est passés à une autre espèce donc on écrit la POcc pour la précedente
				s = Ka(line[3],a,b)*TF/(1+Ka(line[3],a,b)*TF) #on recommence à calculer la POcc pour la nouvelle espèce
				current_gene = line[0] 
		pocc_scores.write(str(s).encode('utf-8')+"\n") 
		
	print('POcc computation done')



# Arguments

if __name__ == "__main__": 
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--scores', help = 'scores files produced by ', required = True)
#	parser.add_argument('-s', '--sequences', help = 'path to the fasta file of the sequences to analyze', required=True) 
#	parser.add_argument('-m', '--matrix', help = 'matrix for sequence analysis', required=True)
	parser.add_argument('-o', '--output', help = 'output path', required = True)
#	parser.add_argument('-n', '--name', help = 'generic name for the different results', required = True)
#	parser.add_argument('-th', '--threshold', help = 'threshold for score computation')
	parser.add_argument('-a', '--a', help = '', default=1.0)
	parser.add_argument('-b', '--b', help = '', default=0.0)
	args = parser.parse_args()
	main(args.scores, args.output,args.a,args.b)



