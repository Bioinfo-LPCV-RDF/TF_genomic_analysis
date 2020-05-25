# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
import scipy.stats as stats
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import argparse


matplotlib.rc('font', serif='Arial') 
matplotlib.rcParams.update({'font.size': 6})

gene_tot=28495


parser = argparse.ArgumentParser()                                               
parser.add_argument("--DAP", "-DAP" , type=int)
parser.add_argument("--ChIP", "-ChIP" , type=int)
parser.add_argument("--common", "-common" , type=int)
parser.add_argument("--name_DAP", "-name_DAP" , type=str)
parser.add_argument("--name_ChIP", "-name_ChIP" , type=str)
parser.add_argument("--out", "-out" , type=str)
parser.add_argument("--col_seq", "-col" , nargs='+',type = int)


args = parser.parse_args()

name_DAP = args.name_DAP
name_ChIP = args.name_ChIP
DAP = args.DAP
ChIP = args.ChIP
common = args.common
out = args.out
col_seq = args.col_seq

# bleu, orange, rouge, vert, gris
palette=("#40a5c7","#f0875a","#f9626e","#08af8e","#bdc3c6","#0D7294","#BD5427","#FF8D99","#08AF8E")


pal_sel=(palette[col_seq[0]],palette[col_seq[1]])

p=stats.hypergeom.sf(common-1,gene_tot,DAP+common,ChIP+common)
p=-(int(np.floor(np.log10(p))))


plt.figure()
plt.figure(figsize=(2.2,2.2))
v=venn2(subsets = (ChIP, DAP, common), set_labels = (name_ChIP, name_DAP),alpha = 1)
v.get_patch_by_id('10').set_color(pal_sel[0])
v.get_patch_by_id('11').set_color(palette[4])
v.get_patch_by_id('01').set_color(pal_sel[1])
v.get_label_by_id('10').set_color("white")
v.get_label_by_id('01').set_color("white")


v.get_label_by_id('01').set_weight("bold")
v.get_label_by_id('10').set_weight("bold")
v.get_label_by_id('11').set_weight("bold")

c=venn2_circles(subsets = (ChIP, DAP, common), linewidth=1.5, color="white")
plt.title('-log(p) = '+str(p))

plt.gca().set_axis_off()
plt.margins(0,0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.savefig(out,bbox_inches = "tight",pad_inches=0)

