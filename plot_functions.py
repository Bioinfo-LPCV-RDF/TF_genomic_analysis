# -*- coding: utf-8 -*-                                                                                                                                                                                            


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import numpy as np
import sys
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
import math
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import re

def divide(a, b):
    if b == 0:
        return np.nan
    else: 
        return a/b

def negative_sets_histo (Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command) :
	fig = plt.figure(1,figsize= (18,10))
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax1 = fig.add_subplot(1,4,1)
	ax1.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax1.axis([0, 3, 0, 2])
	ax1.bar(indexes1, rate , width , color = 'green')
	ax1.set_ylabel('$\sum DR$n_+, ER$n_+, IR$n_+ / $\sum DR$n_-, ER$n_-, IR$n_- , respectively', color = 'green', size = 16)
	for a,b in zip(relative_DR,relative_DR_neg) :
		ax1 = fig.add_subplot(1,4,2)
		ax1.set_xlabel("base pairs between direct repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('DR$n_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("DR$n_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("DR$n_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "D$Rn_+$ = DRn number in the bound set")
		plt.text(2, 4, "D$Rn_-$ = DRn number in the unbound set")
		plt.text(2, 3, "DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		plt.text(2, 6, sys.argv)
		l = plt.axhline(y = 1)
		ax1.legend()
		ax2.legend()
		
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		ax1 = fig.add_subplot(1,4,3)
		ax1.set_xlabel("base pairs between everted repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('E$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("E$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("E$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "ER$n_+$ = ERn number in the bound set")
		plt.text(2, 4, "ER$n_-$ = ERn number in the unbound set")
		plt.text(2, 3, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		ax1.legend()
		ax2.legend()
	for a,b in zip(relative_IR,relative_IR_neg) :	
		ax1 = fig.add_subplot(1,4,4)
		ax1.set_xlabel("base pairs between inverted repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('I$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("I$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("I$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "IR$n_+$ = DRn number in the bound set")
		plt.text(2, 4, "IR$n_-$ = DRn number in the unbound set")
		plt.text(2, 3, "IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		ax1.legend()
		ax2.legend()
	
	plt.show()
		
def negative_sets_curve(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command) :
	fig = plt.figure(1,figsize= (18,10))
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax1 = fig.add_subplot(1,4,1)
	ax1.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax1.axis([0, 3, 0, 2])
	ax1.bar(indexes1, rate , width , color = 'green')
	ax1.set_ylabel('$\sum DR$n_+, ER$n_+, IR$n_+ / $\sum DR$n_-, ER$n_-, IR$n_- , respectively', color = 'green', size = 16)
	plt.text(2, 3, command)
	for a,b in zip(relative_DR,relative_DR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,2)
		ax1.set_xlabel("base pairs between direct repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2,label=r"")
		ax1.set_ylabel("DR$n_+ frequence / DR$n_- frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend()
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,3)
		ax1.set_xlabel("base pairs between everted repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2, label= r"threshold : "+str(c))
		ax1.set_ylabel("$ERn_+$ frequence / $ERn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		plt.legend()
	for a,b in zip(relative_IR,relative_IR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,4)
		ax1.set_xlabel("base pairs between inverted repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2,label=r"")
		ax1.set_ylabel("$IRn_+$ frequence / $IRn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend()
	plt.legend()
	plt.show()
	
def negative_sets_points(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command,output,Interdistance_minValue,load=False,log=False) :	
        matplotlib.rcParams['font.sans-serif']=['Arial']
        matplotlib.rcParams['font.family']="sans-serif"
        SMALL_SIZE=5.4
        # if log==True:
        #     plt.rc('text', usetex = True)
        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
	fig = plt.figure(1,figsize= (3.54,2.2))
	# fig.suptitle(command, fontsize = 11, fontweight='bold')
	gs = gridspec.GridSpec(1, 4,wspace = 0.25, width_ratios=[1, 2,2,2],left=0.1,right=0.97) 
	indexes1 = range(Interdistance_minValue,Interdistance_maxValue + 1)
	width = 1
	ax0 = plt.subplot(gs[0])
        ax0.spines["top"].set_linewidth(0.3)
        ax0.spines["bottom"].set_linewidth(0.3)
        ax0.spines["right"].set_linewidth(0.3)
        ax0.spines["left"].set_linewidth(0.3)
        ax0.tick_params(length=1.5,width=0.5,pad=1.5)
	ax0.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 , ('DR', 'ER', 'IR'))
	turn = 0
        if not load:
            points = ['lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black']
            ax0.set_color_cycle(['lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            print "salut"
            points = ['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue']
            ax0.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
	maxi = []
	for a in rate :
		maxi.append(max(a))
	m = max(maxi)
        if log == False:
			ax0.axis([-1, 3 , 0, m + 0.5])
        else :
			ax0.axis([-1, 3 , -(np.log10(m) + 0.5), np.log10(m) + 0.5])
	for a in rate :
            if log==False :
                plt.plot(indexes1, a, points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
            else :
                plt.plot(indexes1, np.log10(a), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
	    turn = turn + 1
	ax0.set_ylabel("Absolute enrichment", fontsize = 5.4)
	turn = 0
	#plt.text(-3.6, -0.2, command)	
	maxi = []
	for a,b in zip(relative_DR,relative_DR_neg) :
		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_ER,relative_ER_neg) :
		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_IR,relative_IR_neg) :
		maxi.append(max(map(divide, a, b)))
	for n,i in enumerate(maxi): # add by JL, in case too many "nan" in the maxi list
		if math.isnan(i):
			maxi[n] = 0.0
	m = max(maxi)		
	print(m)
	ax1 = plt.subplot(gs[1])
        ax1.spines["top"].set_linewidth(0.3)
        ax1.spines["bottom"].set_linewidth(0.3)
        ax1.spines["right"].set_linewidth(0.3)
        ax1.spines["left"].set_linewidth(0.3)
        ax1.tick_params(length=1.5,width=0.5,pad=1.5)
        ax1.yaxis.set_visible(False)
        if not load:
            ax1.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            ax1.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
	print(len(relative_DR))
	for a,b in zip(relative_DR,relative_DR_neg) :
		indexes1 = np.arange(Interdistance_minValue,Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(2,2,2)
		ax1.set_xlabel("bp between DRs", fontsize = 5.4)
                if log == False:
                    ax1.axis([-1+Interdistance_minValue, Interdistance_maxValue + 1 , 0, m + 0.5])
                else :
                    ax1.axis([-1+Interdistance_minValue, Interdistance_maxValue + 1 , -(np.log10(m) + 0.5), np.log10(m) + 0.5])
		#plt.Circle((indexes1, map(divide, a, b)), color = 'g', radius=50, fill=True)
                if log==False :
                    plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
                else :
                    plt.plot(indexes1, np.log10(map(divide, a, b)), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
		    #plt.title('DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
		
		# t = plt.text(0.5,0.5, 'DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$',horizontalalignment='center', style='italic',verticalalignment='center', transform = ax1.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.text(2, 3.5, "DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
                if log == False :
                    l = plt.axhline(y = 1,linewidth=0.1)
                else :
                    l = plt.axhline(y = 0,linewidth=0.1)
		for i in range (0,Interdistance_maxValue + 1) :
                    if i%5 == 0 :
			plt.axvline(x = i, linewidth=0.3)
		turn = turn + 1
	turn = 0
	ax2 = plt.subplot(gs[2])
        ax2.spines["top"].set_linewidth(0.3)
        ax2.spines["bottom"].set_linewidth(0.3)
        ax2.spines["right"].set_linewidth(0.3)
        ax2.spines["left"].set_linewidth(0.3)
        ax2.tick_params(length=1.5,width=0.5,pad=1.5)        
        ax2.set_ylabel("Normalized enrichment", fontsize = 5.4, labelpad=60)
        if not load:
            ax2.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            ax2.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])

	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	    
		indexes1 = range(Interdistance_minValue,Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(1,7,4)
		ax2.set_xlabel("bp between ERs", fontsize = 5.4)
                if log == False:
                    ax2.axis([-1+Interdistance_minValue, Interdistance_maxValue + 1 , 0, m + 0.5])
                else :
                    ax2.axis([-1+Interdistance_minValue, Interdistance_maxValue + 1 , -(np.log10(m) + 0.5), np.log10(m) + 0.5])
                try :
                    int(c)
                    if log==False :
                        plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, label = r"threshold : " + str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)
                    else :
                        plt.plot(indexes1, np.log10(map(divide, a, b)), points[turn], marker='o',linewidth=0.5, label = r"threshold : " + str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)
                except:
                    plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, label = str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)
		#plt.plot(indexes1, map(divide, a, b), points[turn], marker='o', alpha = 0.70, markersize = 13,markeredgewidth = 0.0)
		#plt.figlegend( d,('lightsalmon', 'tomato', 'r'), loc = 'lower center', ncol=5, labelspacing=0. )

		# ax2.set_ylabel("ERn+ frequence / ERn- frequence", size = 16)
		# t = plt.text(0.5,0.5, 'ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', horizontalalignment='center', style='italic',verticalalignment='center', transform = ax2.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.title('ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
		#plt.text(0, 3.5, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
                if log == False :
                    l = plt.axhline(y = 1,linewidth=0.1)
                else :
                    l = plt.axhline(y = 0,linewidth=0.1)                 
		plt.legend(bbox_to_anchor=(0.5,0.5), ncol=(len(threshold)-1//2), mode = "expand", fontsize=5.4 ,borderaxespad = -15)
		turn = turn + 1
		for i in range (0,Interdistance_maxValue + 1) :
                    if i%5 == 0 :
			plt.axvline(x = i, linewidth=0.3)
	turn = 0
	ax3 = plt.subplot(gs[3])
        ax3.spines["top"].set_linewidth(0.3)
        ax3.spines["bottom"].set_linewidth(0.3)
        ax3.spines["right"].set_linewidth(0.3)
        ax3.spines["left"].set_linewidth(0.3)
        ax3.tick_params(length=1.5,width=0.5,pad=1.5)
        if not load:
            ax3.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            ax3.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
	for a,b in zip(relative_IR,relative_IR_neg) :
		indexes1 = range(Interdistance_minValue, Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(1,7,6)
		ax3.set_xlabel("bp between IRs", fontsize = 5.4)
                if log == False:
                    ax3.axis([-1+Interdistance_minValue, Interdistance_maxValue + 1 , 0, m + 0.5])
                else :
                    ax3.axis([-1+Interdistance_minValue, Interdistance_maxValue + 1 , -(np.log10(m) + 0.5), np.log10(m) + 0.5])
                if log==False :
                    plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
                else :
                    plt.plot(indexes1, np.log10(map(divide, a, b)), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
		# ax3.set_ylabel("IRn+ frequence / IRn- frequence", size = 16)
		# t = plt.text(0.5,0.5, 'IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', horizontalalignment='center', style='italic',verticalalignment='center', transform = ax3.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.title('IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
                if log == False :
                    l = plt.axhline(y = 1,linewidth=0.1)
                else :
                    l = plt.axhline(y = 0,linewidth=0.1)
		turn = turn + 1
		for i in range (0,Interdistance_maxValue + 1):
                    if i%5 == 0 :
			plt.axvline(x = i, linewidth=0.3)
	plt.subplots_adjust(top=0.85)
#        plt.tight_layout(fig)
	fig.savefig(output)
        if log==True:
            labels_ax0 = [item.get_text() for item in ax0.get_yticklabels()]
            labels_ax2 = [item.get_text() for item in ax2.get_yticklabels()]
            new_labels_ax0=[]
            new_labels_ax2=[]
            for elt in labels_ax0 :
                # elt=elt.replace('$','R',1)
                # elt=elt.replace('$','}$')
                # elt=elt.replace('R','$10^{')
                # new_labels_ax0.append(elt)
                elt=10**float(elt.encode('utf-8').replace('−','-')) ### use if you do not like scientific expression. 
                new_labels_ax0.append(str(round(elt, 2))) ### You need to remove plt.rc('text', usetex = True) at the beginning of the function
            for elt in labels_ax2 :
                # elt=elt.replace('$','R',1)
                # elt=elt.replace('$','}$')
                # elt=elt.replace('R','$10^{')
                # new_labels_ax2.append(elt)
                elt=10**float(elt.encode('utf-8').replace('−','-')) ### use if you do not like scientific expression. 
                new_labels_ax2.append(str(round(elt, 2))) ### You need to remove plt.rc('text', usetex = True) at the beginning of the function
            ax0.set_yticklabels(new_labels_ax0)
            ax2.set_yticklabels(new_labels_ax2)
            # ax0.yaxis.set_tick_params(pad=-2)
            # ax2.yaxis.set_tick_params(pad=-2)
            ax0.set_ylabel("Absolute enrichment (log)", fontsize = 5.4,labelpad=-1)
            ax2.set_ylabel("Normalized enrichment (log)", fontsize = 5.4,labelpad= 59)
            fig.savefig(output)

        
def negative_sets_points_ER(Interdistance_maxValue,relative_R,relative_R_neg,threshold,rate,command,output,Interdistance_minValue,load=False,log=False,no_threshold_panel=False,no_absolute_panel=False,maxy=0) :	
        matplotlib.rcParams['font.sans-serif']=['Arial']
        matplotlib.rcParams['font.family']="sans-serif"
        SMALL_SIZE=5.4
        # if log==True:
        #     plt.rc('text', usetex = True)
        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
	fig = plt.figure(1,figsize= (3.54,2.2))
	# fig.suptitle(command, fontsize = 11, fontweight='bold')
        if not load:
            points = ['#40A5C7', '#307C95', '#205364', '#102932', 'black'] #['lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black']
        else:
            points = ['#40A5C7', '#307C95', '#205364', '#102932', 'black'] #['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue']
        if no_absolute_panel is False :
	    gs = gridspec.GridSpec(1, 2,wspace = 0.25, width_ratios=[1, 12],left=0.1,right=0.97) 
	    indexes1 = range(Interdistance_maxValue + 1)
	    width = 1
	    ax0 = plt.subplot(gs[0])
            if not load:
                #ax0.set_color_cycle(['lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
				ax0.set_color_cycle([ '#40A5C7', '#307C95', '#205364', '#102932', 'black'])
            else:
                #ax0.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
				ax0.set_color_cycle([ '#40A5C7', '#307C95', '#205364', '#102932', 'black'])
            ax0.spines["top"].set_linewidth(0.3)
            ax0.spines["bottom"].set_linewidth(0.3)
            ax0.spines["right"].set_linewidth(0.3)
            ax0.spines["left"].set_linewidth(0.3)
            ax0.tick_params(length=1.5,width=0.5,pad=1.5)
	    ax0.set_xlabel("", size = 16)
	    indexes1 = np.arange(1)
	    plt.xticks(indexes1 , ('R',))
	    turn = 0
            rateER=[]
            for a in rate :
                rateER.append(a[1])
            rate=list(rateER)
	    m = max(rate)
            if log == False:
				ax0.axis([-2, 2 , 0, m + 0.5])
            else :
				print("lol2")
				ax0.axis([-2, 2 , -(np.log10(m) + 0.2), np.log10(m) + 0.2])
            for a in rate :
                if log==False :
                    plt.plot(indexes1, a, points[turn], marker='o',linewidth=1.25, alpha = 0.7, markersize = 3, markeredgewidth = 0.0)
                else :
                    plt.plot(indexes1, np.log10(a), points[turn], marker='o',linewidth=1.25, alpha = 0.7, markersize = 3, markeredgewidth = 0.0)
                turn+=1
	    ax0.set_ylabel("Absolute enrichment", fontsize = 5.4)
            gs_ind=0
            ##  normalized enrichment panel form here
        else :
            gs = gridspec.GridSpec(1, 1,left=0.1,right=0.97)
            gs_ind=1            
	turn = 0
	#plt.text(-3.6, -0.2, command)	
	maxi = []
#	for a,b in zip(relative_DR,relative_DR_neg) :
#		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_R,relative_R_neg) :
		maxi.append(max(map(divide, a, b)))
#	for a,b in zip(relative_IR,relative_IR_neg) :
#		maxi.append(max(map(divide, a, b)))
	for n,i in enumerate(maxi): # add by JL, in case too many "nan" in the maxi list
		if math.isnan(i):
			maxi[n] = 0.0
	print(maxi)
	m = max(maxi)		
	print(m)
	turn = 0
	ax2 = plt.subplot(gs[1-gs_ind])
        ax2.spines["top"].set_linewidth(0.3)
        ax2.spines["bottom"].set_linewidth(0.3)
        ax2.spines["right"].set_linewidth(0.3)
        ax2.spines["left"].set_linewidth(0.3)
        ax2.tick_params(length=1.5,width=0.5,pad=1.5)        
        ax2.set_ylabel("Normalized Enrichment", fontsize = 8)
        if not load:
            #ax2.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
			ax2.set_color_cycle([ '#40A5C7', '#307C95', '#205364', '#102932', 'black'])
			# '#F0875A', '#B46544', '#78442D' SEP3
			# '#40A5C7', '#307C95', '#205364' SEP3AG
			# '#F9626E', '#BB4A52', '#7C3137' SEP3delAG
        else:
            #ax2.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
			ax2.set_color_cycle([ '#40A5C7', '#307C95', '#205364', '#102932', 'black'])
	for a, b,c in zip(relative_R,relative_R_neg,threshold)  :
                # print(map(divide, a, b))                        
                indexes1 = range(Interdistance_minValue,Interdistance_maxValue + 1) #+ Interdistance_minValue
		#ax1 = fig.add_subplot(1,7,4)
		ax2.set_xlabel("bp between CArG boxes", fontsize = 8)
                if log == False:
					ylimit = max(m + 0.5, 3.0)
					print(ylimit)
					print("lol")
					ax2.axis([-1 + Interdistance_minValue, Interdistance_maxValue + 1 , 0, ylimit])
					if maxy != 0:
						ax2.axis([-1 + Interdistance_minValue, Interdistance_maxValue + 1 , 0, maxy])
                else :
                    ax2.axis([-1 + Interdistance_minValue, Interdistance_maxValue + 1 , -(np.log10(m) + 0.5), np.log10(m) + 0.5])
                    if maxy != 0:
                        ax2.axis([-1 + Interdistance_minValue, Interdistance_maxValue + 1 , -np.log10(maxy), np.log10(maxy)])
                    # ax2.get_yaxis().set_major_formatter(FuncFormatter(lambda x,pos :"10**"x))
                    # ax2.axis([-1, Interdistance_maxValue + 1 , 0, m + 0.5])
                    # ax2.set_yscale('log')                    
                try :
                    int(c)
                    if log==False :
                        plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=1.25, label = r"threshold : " + str(c), alpha = 0.7, markersize = 3,markeredgewidth = 0.0)
                    else :
                        plt.plot(indexes1, np.log10(map(divide, a, b)), points[turn], marker='o',linewidth=1.25, label = r"threshold : " + str(c), alpha = 0.7, markersize = 3,markeredgewidth = 0.0)
                        # plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, label = r"threshold : " + str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)

                except:
                    plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, label = str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)
                    
		#plt.plot(indexes1, map(divide, a, b), points[turn], marker='o', alpha = 0.70, markersize = 13,markeredgewidth = 0.0)
		#plt.figlegend( d,('lightsalmon', 'tomato', 'r'), loc = 'lower center', ncol=5, labelspacing=0. )

		# ax2.set_ylabel("ERn+ frequence / ERn- frequence", size = 16)
		# t = plt.text(0.5,0.5, 'ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', horizontalalignment='center', style='italic',verticalalignment='center', transform = ax2.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.title('ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
		#plt.text(0, 3.5, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
                if log == False :
                    l = plt.axhline(y = 1,linewidth=0.1, color="black", linestyle=":")
                else :
                    l = plt.axhline(y = 0,linewidth=0.1, color="black", linestyle=":")
                if no_threshold_panel == False:
		    plt.legend(bbox_to_anchor=(0.5,0.5), ncol=(len(threshold)-1//2), mode = "expand", fontsize=5.4 ,borderaxespad = -15)
		turn = turn + 1
#		for i in range (0,Interdistance_maxValue + 1) :
#                    if i%10 == 0 :
#			plt.axvline(x = i, linewidth=0.3,alpha=0.1)
	plt.subplots_adjust(top=0.85)
	fig.savefig(output)
        if log==True:
            if no_absolute_panel is False :
                labels_ax0 = [item.get_text() for item in ax0.get_yticklabels()]
                new_labels_ax0=[]
                for elt in labels_ax0 :
                    # elt=elt.replace('$','R',1)
                    # elt=elt.replace('$','}$')
                    # elt=elt.replace('R','$10^{')
                    # new_labels_ax0.append(elt)
                    elt=10**float(elt.encode('utf-8').replace('−','-')) ### use if you do not like scientific expression. 
                    new_labels_ax0.append(str(round(elt, 2))) ### You need to remove plt.rc('text', usetex = True) at the beginning of the function
                ax0.set_yticklabels(new_labels_ax0)
                ax0.set_ylabel("Absolute enrichment (log)", fontsize = 5.4,labelpad=-1)
            labels_ax2 = [item.get_text() for item in ax2.get_yticklabels()]
            new_labels_ax2=[]
            for elt in labels_ax2 :
                # elt=elt.replace('$','R',1)
                # elt=elt.replace('$','}$')
                # elt=elt.replace('R','$10^{')
                # new_labels_ax2.append(elt)
                elt=10**float(elt.encode('utf-8').replace('−','-')) ### use if you do not like scientific expression. 
                new_labels_ax2.append(str(round(elt, 2))) ### You need to remove plt.rc('text', usetex = True) at the beginning of the function
            ax2.set_yticklabels(new_labels_ax2)
            # ax0.yaxis.set_tick_params(pad=-2)
            # ax2.yaxis.set_tick_params(pad=-2)
            ax2.set_ylabel("Normalized enrichment (log)", fontsize = 8,labelpad= 2)
            fig.savefig(output)


