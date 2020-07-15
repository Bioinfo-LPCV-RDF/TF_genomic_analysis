# -*- coding: utf-8 -*-
#!usr/bin/python
################################################                General intels                ################################################                
#description     : lui demander (ou à Jérémy)
#author          :Arnaud Stigliani
#date            : 01/11/2018
#version         :0.1
#usage           : [script name] -h
#python_version  :2.7.9
################################################                import packages                ################################################    
# library
import os
import metaseq
import multiprocessing
from pybedtools import BedTool
import numpy as np
from matplotlib import pyplot as plt
from math import log10
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import gridspec
import argparse
processes = multiprocessing.cpu_count()




################################################                Functions                ################################################
#def read_table(table,column_to_get):
#	table_col=[]
#	with open(table) as IN:
#		IN.next()
#		for line in IN:
#			columns = line.split("\r")[0].split("\n")[0].split("\t")
#			value = float(columns[column_to_get])
#			if value > 0.0:
#				table_col.append(log10(value))
#			else:
#				table_col.append(log10(0.1))
#	return np.array(table_col)
################################################                MAIN                ################################################
def main(data_dir, peak_file, results, sortBy, list_name):

	regions=np.loadtxt(peak_file,dtype={'names': ('chr','start','stop','file'),'formats':('S15','int','int','S15')})['file']
	peaks=BedTool(peak_file)
#	if not os.path.exists(results+"_".join(list_name)+'.npz'):
	print "generating"
	ip_name1 = peak_file.split("/")[-1].split("_")[0]+".bw"
	ip_name2 = peak_file.split("/")[-1].split("_")[1].split("_processed")[0]+".bw"
	print(ip_name1,ip_name2)
	ip_signal_1 = metaseq.genomic_signal(os.path.join(data_dir,ip_name1),'bigwig')
	ip_signal_2 = metaseq.genomic_signal(os.path.join(data_dir,ip_name2),'bigwig')
	ip_array_1 = ip_signal_1.array( peaks, bins=800,processes=processes,method="get_as_array")
	ip_array_2 = ip_signal_2.array( peaks, bins=800,processes=processes,method="get_as_array")
	metaseq.persistence.save_features_and_arrays(features=peaks,arrays={'ip1': ip_array_1,'ip2': ip_array_2},prefix=results+"_".join(list_name),link_features=True,overwrite=True)

	features, arrays = metaseq.persistence.load_features_and_arrays(prefix=results+"_".join(list_name))
	ip1_mean=np.apply_along_axis(np.sum,1,arrays['ip1'][::,300:500])
	ip2_mean=np.apply_along_axis(np.sum,1,arrays['ip2'][::,300:500])
	#print(arrays['ip2'])

	#CFR=(ip1_mean)/(ip2_mean)
	#print(CFR)
	if sortBy == 1:
		ip1_mean_order=np.argsort(ip1_mean)
		path_fig = results+"heatmaps_"+list_name[0]+"_"+list_name[1]+"_sorted_on_"+list_name[0]+".png"
	elif sortBy == 2:
		ip1_mean_order=np.argsort(ip2_mean)
		path_fig = results+"heatmaps_"+list_name[0]+"_"+list_name[1]+"_sorted_on_"+list_name[1]+".png"
	#elif sortBy == 3:
		#ip1_mean_order=np.argsort(CFR)
		#path_fig = results+"heatmaps_"+list_name[0]+"_"+list_name[1]+"_sorted_on_CFR.png"
	
	#print(ip1_mean_order)
	ip1=arrays['ip1'][ip1_mean_order,:][::1]
	ip2=arrays['ip2'][ip1_mean_order,:][::1]
	#CFR_ordered=CFR[ip1_mean_order][::1]
	
	x = np.linspace(-400, 400, 800)
	regions_ordered=regions[ip1_mean_order][::1]
	Vmax=max([ip1.mean(),ip2.mean()])*2
	Vmin=min([ip1.min(),ip2.min()])
	gs = gridspec.GridSpec(2, 2, width_ratios=[8,8],height_ratios=([8,2]))
	gs.update(wspace=0.02, hspace=0.02) 
	
	plt.rcParams['font.family'] = 'Arial'
	plt.rcParams['font.size'] = 10
	fig=plt.figure(figsize=(9,5))
	ax1=fig.add_subplot(gs[0,0])
	ax1.pcolormesh(ip1,vmin=Vmin,vmax=Vmax,cmap=cm.Reds)
	ax1.set_title(list_name[0]+' coverage')
	ax1.xaxis.set_ticklabels([])
	ax1.xaxis.set_ticks_position('none') 
	ax2=fig.add_subplot(gs[0,1])
	ax2.pcolormesh(ip2,vmin=Vmin,vmax=Vmax,cmap=cm.Reds)
	ax2.set_title(list_name[1]+' coverage')
	ax2.yaxis.tick_right()
	ax2.xaxis.set_ticklabels([])
	ax2.xaxis.set_ticks_position('none') 
	#ax3=fig.add_subplot(gs[2])
##	CFR_min=CFR_interval[0]
##	CFR_max=CFR_interval[1]
	
	#CFR_min=CFR.min()
	#CFR_max=CFR.max()
	#pcm=ax3.pcolormesh(np.column_stack(((CFR_ordered),(CFR_ordered))),vmin=CFR_min,vmax=3,cmap=cm.Blues)
	#ax3.set_title('CFR (W='+str(round(CFR_min,2))+', B='+str(round(CFR_max,2))+')')
	#ax3.yaxis.set_ticklabels([])
	#ax3.xaxis.set_ticklabels([])
	
	
#	ax4=fig.add_subplot(gs[3])
#	V_min=Dnase_leaf_ordered.min()
#	V_max=Dnase_leaf_ordered.max()
#	pcm=ax4.pcolormesh(np.column_stack(((Dnase_leaf_ordered),(Dnase_leaf_ordered))),vmin=V_min,vmax=V_max,cmap=cm.Blues)
#	ax4.set_xlabel('Dnase Leaf (W='+str(round(V_min,2))+', B='+str(round(V_max,2))+')')
#	ax4.yaxis.set_ticklabels([])
#	ax4.xaxis.set_ticklabels([])
	
	ax3=fig.add_subplot(gs[1,0])
	ax3.plot(x,ip1.mean(axis=0))
	
	ax4=fig.add_subplot(gs[1,1])
	ax4.plot(x,ip2.mean(axis=0))
	ax4.yaxis.tick_right()
	#fig.colorbar(pcm, extend='max')
	#plt.tight_layout()






	plt.savefig(path_fig)







## Create a figure and axes
#fig = plt.figure()
#ax = fig.add_subplot(111)


## Plot the IP:
#ax.plot(
    ## use the x-axis values we created
    #x,

    ## axis=0 takes the column-wise mean, so with
    ## 100 columns we'll have 100 means to plot
    #arrays['ip'].mean(axis=0),

    ## Make it red
    #color='r',

    ## Label to show up in legend
    #label='IP')


## Do the same thing with the input
#ax.plot(
    #x,
    #arrays['input'].mean(axis=0),
    #color='k',
    #label='input')


## Add a vertical line at the TSS, at position 0
#ax.axvline(0, linestyle=':', color='k')


## Add labels and legend
#ax.set_xlabel('Distance from TSS (bp)')
#ax.set_ylabel('Average read coverage (per million mapped reads)')
#ax.legend(loc='best');







	########## debug #######
        ########################

#	fig=plt.figure(figsize=(9,5))
#	regions[np.where(regions=='35S')]="blue"
#	regions[np.where(regions=='LFY-FL')]="red"
#	regions[np.where(regions=='both')]="green"
#	plt.scatter(ip1_mean,ip2_mean,s=5,alpha=0.2,c=regions)
#	plt.xlim(0.1,10000)
#	plt.ylim(0.1,10000)
#	plt.xscale('log')
#	plt.yscale('log')
#	plt.tight_layout()
#	plt.savefig(path_fig+"lol.png")
	########################
	########################

################################################                Arguments                ################################################ 
if __name__ == "__main__": 
	parser = argparse.ArgumentParser() 

	parser.add_argument("--bamdir","-b",required =True)
	parser.add_argument("--peakfile","-p",required =True)
	parser.add_argument("--resultdir","-r",required =True)
#	parser.add_argument("--CFR","-c",default =[1,2],nargs="+",type=int)
	parser.add_argument("--sort","-s", type=int, choices=[1,2,3,4],help = "1: sort on TF1; 2: sort on TF2; 3: sort on CFR",required =True)
	parser.add_argument("--name","-n", nargs=2, help="name of the samples provided, in same order as in the peak file",required =True)
#	parser.add_argument("--all","-a",action='store_true', default= False)
	args = parser.parse_args()
	data_dir = args.bamdir
	peak_file = args.peakfile
	results = args.resultdir
	sortBy = args.sort
	list_name = args.name
	main(data_dir, peak_file, results, sortBy, list_name)




#data_dir="../results/mapping/DAP_qsub_dovetail/"
#peak_file="../results/analyses/LIB5-1_001_LIB5-3_001_processed_for_heatmap.bed"
#results="../results/analyses/"


