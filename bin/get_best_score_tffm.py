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

parser.add_argument("--output", "-o", type=str)
parser.add_argument("--fasta_pos", "-pos", type=str)
parser.add_argument("--tffm", "-t", type=str)

args = parser.parse_args()
output = args.output
tffm=args.tffm
fasta_pos=args.fasta_pos


hit_pos_list=[]
tffm_first_order = tffm_module.tffm_from_xml(tffm, TFFM_KIND.FIRST_ORDER)

print tffm
print fasta_pos

for hit_pos in tffm_first_order.scan_sequences(fasta_pos, only_best=True):
    hit_pos_list.append(str(hit_pos))
    
    
with open(output,"w") as f1:
    for hit_pos in hit_pos_list:
        f1.write(hit_pos+"\n")        
