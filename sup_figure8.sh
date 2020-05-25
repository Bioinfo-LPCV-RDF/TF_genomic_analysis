#!/bin/bash

data=../data
dirout=../results
mkdir -p $dirout
tair10=$data/tair10.fas

###### ChIP data #########

ChIP_SEP3=$data/ChIP_SEP3.bed
ChIP_AG=$data/ChIP_AG.bed

ChIP_SEP3AG=$dirout/ChIP_SEP3AG.bed
#intersect ChIP SEP3 x ChIP AG (peaks)
bedtools intersect -a $ChIP_SEP3 -b $ChIP_AG -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $1,$2,$3}' | uniq > $ChIP_SEP3AG

###### DAP data #########

DAP_SEP3AG=$data/DAP_SEP3AG.bed

#### intersect ChIP DAP ####
ChIP_only=$dirout/ChIP_only.bed
DAP_only=$dirout/DAP_only.bed
common=$dirout/common.bed

bedtools intersect -a $ChIP_SEP3AG -b $DAP_SEP3AG > $common

bedtools intersect -v -a $ChIP_SEP3AG -b $common > $ChIP_only
bedtools intersect -v -a $DAP_SEP3AG  -b $common > $DAP_only

DAP_only_fas=$dirout/DAP_only.fas
bedtools getfasta -fi $tair10 -fo  $DAP_only_fas -bed $DAP_only

ChIP_only_fas=$dirout/ChIP_only.fas
bedtools getfasta -fi $tair10 -fo  $ChIP_only_fas -bed $ChIP_only


#### get scores ########
SEP3AG_tffm=$data/tffm_first_order.xml
bs_ChIP=$dirout/bs_ChIP
bs_DAP=$dirout/bs_DAP

python get_best_score_tffm.py -pos $ChIP_only_fas -o $bs_ChIP -t $SEP3AG_tffm
python get_best_score_tffm.py -pos $DAP_only_fas -o $bs_DAP -t $SEP3AG_tffm

Rscript plot_scores_boxplots.r $bs_ChIP $bs_DAP
