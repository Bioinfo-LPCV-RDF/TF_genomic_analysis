# coding: utf-8

from os import system
from os import mkdir
from os.path import isdir
from os.path import basename
from os.path import dirname
from os.path import join
import sys
sys.path.append("/home/304.3-STRUCTPDEV/lib/TFFM")
import tffm_module
from constants import TFFM_KIND
import argparse

parser = argparse.ArgumentParser()                                               

parser.add_argument("--results", "-r", type=str)
parser.add_argument("--fasta", "-f", type=str)
parser.add_argument("--matrix", "-m", type=str)
parser.add_argument("--place", "-p", type=int)
args = parser.parse_args()
results = args.results
meme_motif=args.matrix
fasta=args.fasta
place=args.place
print fasta
print results
print meme_motif
print place

tffm_first_order = tffm_module.tffm_from_meme(meme_motif, TFFM_KIND.FIRST_ORDER,place)
out = open(results+"/tffm_summary_initial.svg", "w")
tffm_first_order.print_summary_logo(out)
out.close()
tffm_first_order.write(results+"/tffm_first_order_initial.xml")
out = open(results+"/tffm_dense_initial.svg", "w")
tffm_first_order.print_dense_logo(out)
out.close()
tffm_first_order.train(fasta)
tffm_first_order.write(results+"/tffm_first_order.xml")
out = open(results+"/tffm_summary.svg", "w")
tffm_first_order.print_summary_logo(out)
out.close()
out = open(results+"/tffm_dense.svg", "w")
tffm_first_order.print_dense_logo(out)
out.close()









