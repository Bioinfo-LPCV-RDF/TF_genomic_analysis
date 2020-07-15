# -*- coding: utf-8 -*-

from os.path import basename
from os.path import join
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f1","--file1", type=str, help='file 1')
parser.add_argument("-f2","--file2", type=str, help='file 2')
parser.add_argument("-o","--output", type=str, help='output directory')
parser.add_argument("-d","--data", type=str, default="", help='data directory')

args = parser.parse_args()

file1 = args.file1
file2 = args.file2
output = args.output
data=args.data


if data == "":
	data = output+"/"+file1+"_"+file2+"/"
else:
	data=data+"/"
	
print(data)


all_peaks=[]
with open(data+file1+"_"+file2+"_peaks.bed","r") as f1:
	for line in f1:
		line=line.strip("\n")
		l=line.split("\t")
		l[1]=int(l[1])
		l[2]=int(l[2])
		l.append(False)
		all_peaks.append(l)
   

	 
i=0
j=0
k=0
while (i < len(all_peaks) -1):
	j=i+1
	
	while (j < len(all_peaks) -1):
		# Start J >= start I
		if all_peaks[i][0] == all_peaks[j][0] and all_peaks[i][4] is not True  : # même CHR
			if  all_peaks[j][2] <= all_peaks[i][2] : # stop de J avant stop I
				all_peaks[j][4]=True
				all_peaks[i][4]=True
				all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3]])
				if (all_peaks[j][1] - all_peaks[i][1]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 : 
					all_peaks.insert(i,[all_peaks[i][0],all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],False])
					i+=1
					j+=1
				if (all_peaks[i][2] - all_peaks[j][2]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 :
					all_peaks.insert(j+1,[all_peaks[i][0],all_peaks[j][2]+1,all_peaks[i][2],all_peaks[i][3],False])
			elif ((all_peaks[j][1] < all_peaks[i][2]) and (all_peaks[j][2] > all_peaks[i][2])) : # start de J contenu dans I
				if all_peaks[i][2] -all_peaks[i][1] < all_peaks[j][2] -all_peaks[j][1]: # J plus petit que I
					if (all_peaks[i][2] - all_peaks[j][1]) /float(all_peaks[i][2] - all_peaks[i][1]) > 0.8:
						all_peaks[j][4]=True
						all_peaks[i][4]=True
						all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3]])
						if (all_peaks[i][2] - all_peaks[j][2]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 :
							all_peaks.insert(j+1,[all_peaks[i][0],all_peaks[j][2]+1,all_peaks[i][2],all_peaks[i][3],False])
				else: # J plus grand ou égal à I
					if (all_peaks[i][2] - all_peaks[j][1]) /float(all_peaks[j][2] - all_peaks[j][1]) > 0.8:
						all_peaks[j][4]=True
						all_peaks[i][4]=True
						all_peaks[i].append([all_peaks[j][0],all_peaks[j][1],all_peaks[j][2],all_peaks[j][3]])
						if (all_peaks[j][1] - all_peaks[i][1]) / float(all_peaks[i][2] - all_peaks[i][1] ) > 0.5 :
							all_peaks.insert(i,[all_peaks[i][0],all_peaks[i][1],all_peaks[j][1]-1,all_peaks[i][3],False])
							i+=1
							j+=1						
			elif ((all_peaks[j][1] >= all_peaks[i][2])): #J après I
				k+=1
				break
			# la
		else : # CHR différents
			break
		# if j-i == 2:
		#	 print j-i,all_peaks[i][0],all_peaks[i][1],all_peaks[i][2],all_peaks[j][1],all_peaks[j][2]
		j+=1
	i+=1
	if  i%1000 ==0 :
		 print(i)
		
print("k="+str(k))
print("i="+str(i))
print("j="+str(j))


with open(data+file1+"_"+file2+"_peaks_uniques.bed","w") as f1:
	for elt in all_peaks:
		if elt[4] is False:
			f1.write(str(elt[0])+"\t"+str(elt[1])+"\t"+str(elt[2])+"\t"+str(elt[3])+"\n")

with open(data+file1+"_"+file2+"_peaks_merged.bed","w") as f1:
	for elt in all_peaks:
		if len(elt) > 5:
			mini=elt[1]
			maxi=elt[2]
			regions=elt[3]
			for elt2 in elt:
				if type(elt2) == list:
					regions="both" # str(regions)+","+str(elt2[3])
					if elt2[1]>mini:
						mini=elt2[1]
						if elt2[2]<maxi:
							maxi=elt2[2]
			f1.write(str(elt[0])+"\t"+str(mini)+"\t"+str(maxi)+"\t"+str(regions)+"\n")
		
