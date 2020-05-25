#!/bin/bash

dirout=../results/figure5
mkdir -p $dirout

###### regulated_gene_lists #########
AG_reg=../data/DEG_AG
SEP3_reg=../data/DEG_SEP3
SEP3AG_over_under=../data/DEG_SEP3_and_AG.tsv

SEP3AG_reg=$dirout/SEP3AG_reg

awk '{print $1}' $SEP3AG_over_under > $SEP3AG_reg


###### ChIP data #########
chip_SEP3_preprocessed=../data/ChIP_SEP3.bed
chip_AG_preprocessed=../data/ChIP_AG.bed

chip_SEP3=$dirout/ChIP_SEP3.bed
chip_AG=$dirout/ChIP_AG.bed

awk -v OFS="\t" '{print $1,$2,$3}' $chip_SEP3_preprocessed > $chip_SEP3
awk -v OFS="\t" '{print $1,$2,$3}' $chip_AG_preprocessed > $chip_AG


chip_SEP3AG=$dirout/Intersect_chipSEP3_x_chipAG.bed

#intersect ChIP SEP3 x ChIP AG (peaks)
bedtools intersect -a $chip_SEP3 -b $chip_AG -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $1,$2,$3}' > $chip_SEP3AG

###### DAP data #########

DAP_SEP3AG=../data/DAP_SEP3AG.bed
DAP_SEP3=../data/DAP_SEP3.bed


###### gene bed file  #########
extendedAT_genes=../data/all_genes_extended.bed

#### get bound genes : ChIP ###
genes_boundSEP3AG_ChIP=$dirout/genes_boundSEP3AG_ChIP
genes_boundSEP3_ChIP=$dirout/genes_boundSEP3_ChIP
genes_boundAG_ChIP=$dirout/genes_boundAG_ChIP

echo "========= ChIP-Seq ========== "

bedtools intersect -a $chip_SEP3AG -b $extendedAT_genes -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $7} ' | sort -u > $genes_boundSEP3AG_ChIP
echo "Genes bound by SEP3/AG "$(wc -l  $genes_boundSEP3AG_ChIP | awk '{print $1}')


bedtools intersect -a $chip_SEP3 -b $extendedAT_genes -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $7} ' | sort -u > $genes_boundSEP3_ChIP
echo "Genes bound by SEP3 "$(wc -l  $genes_boundSEP3_ChIP | awk '{print $1}')

bedtools intersect -a $chip_AG -b $extendedAT_genes -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $7} ' | sort -u > $genes_boundAG_ChIP
echo "Genes bound by AG "$(wc -l  $genes_boundAG_ChIP | awk '{print $1}')


#### get bound genes : DAP ###
genes_boundSEP3AG_DAP=$dirout/genes_boundSEP3AG_DAP
genes_boundSEP3_DAP=$dirout/genes_boundSEP3_DAP

echo "========= DAP-Seq ========== "

bedtools intersect -a $DAP_SEP3AG -b $extendedAT_genes -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $7} ' | sort -u > $genes_boundSEP3AG_DAP
echo "Genes bound by SEP3AG "$(wc -l  $genes_boundSEP3AG_DAP | awk '{print $1}')


bedtools intersect -a $DAP_SEP3 -b $extendedAT_genes -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $7} ' | sort -u > $genes_boundSEP3_DAP
echo "Genes bound by SEP3 "$(wc -l  $genes_boundSEP3_DAP | awk '{print $1}')

declare -A DAP_bound_genes
DAP_bound_genes=( ["SEP3AG"]=$genes_boundSEP3AG_DAP ["SEP3"]=$genes_boundSEP3_DAP)

declare -A color
color=( ["SEP3AG"]="0" ["SEP3"]="1")


for elt in "SEP3AG" "SEP3"
do
    echo "======== DAP-Seq : ${elt}  ==========="
    ##### creating directories #####
    mkdir -p $dirout/DAP_${elt}_vs_ChIP_SEP3AG  $dirout/DAP_${elt}_vs_ChIP_SEP3 $dirout/DAP_${elt}_vs_ChIP_AG $dirout/figures

    # 1 DAP_${elt}_vs_ChIP_SEP3AG
    cat "${DAP_bound_genes[${elt}]}" $genes_boundSEP3AG_ChIP | sort | uniq -c | grep "2\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common
    cat $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common "${DAP_bound_genes[${elt}]}"  | sort | uniq -c | grep "1\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_SEP3AG/DAP
    cat $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common $genes_boundSEP3AG_ChIP  | sort | uniq -c | grep "1\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_SEP3AG/ChIP
    nDAP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG/DAP | awk '{print $1}')
    nChIP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG/ChIP | awk '{print $1}')
    ncommon=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common | awk '{print $1}')
    echo $nDAP $ncommon $nChIP
    python make_venn_figure5.py -DAP $nDAP -ChIP $nChIP -common $ncommon -name_DAP "$(echo -e DAP-Seq\\n${elt})" -name_ChIP "$(echo -e ChIP-Seq\\nSEP3/AG)" -out $dirout/figures/DAP_${elt}_vs_ChIP_SEP3AG.svg -col 5 "${color[${elt}]}"

    # 2 DAP_${elt}_vs_ChIP_SEP3
    cat "${DAP_bound_genes[${elt}]}" $genes_boundSEP3_ChIP | sort | uniq -c | grep "2\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_SEP3/common
    cat $dirout/DAP_${elt}_vs_ChIP_SEP3/common "${DAP_bound_genes[${elt}]}"  | sort | uniq -c | grep "1\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_SEP3/DAP
    cat $dirout/DAP_${elt}_vs_ChIP_SEP3/common $genes_boundSEP3_ChIP  | sort | uniq -c | grep "1\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_SEP3/ChIP
    nDAP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3/DAP | awk '{print $1}')
    nChIP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3/ChIP | awk '{print $1}')
    ncommon=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3/common | awk '{print $1}')
    echo $nDAP $ncommon $nChIP
    python make_venn_figure5.py -DAP $nDAP -ChIP $nChIP -common $ncommon -name_DAP "$(echo -e DAP-Seq\\n${elt})" -name_ChIP "$(echo -e ChIP-Seq\\nSEP3)" -out $dirout/figures/DAP_${elt}_vs_ChIP_SEP3.svg -col 6 "${color[${elt}]}"


    # 3 DAP_${elt}_vs_ChIP_AG
    cat "${DAP_bound_genes[${elt}]}" $genes_boundAG_ChIP | sort | uniq -c | grep "2\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_AG/common
    cat $dirout/DAP_${elt}_vs_ChIP_AG/common "${DAP_bound_genes[${elt}]}"  | sort | uniq -c | grep "1\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_AG/DAP
    cat $dirout/DAP_${elt}_vs_ChIP_AG/common $genes_boundAG_ChIP  | sort | uniq -c | grep "1\sA" | awk '{print $2}' > $dirout/DAP_${elt}_vs_ChIP_AG/ChIP
    nDAP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_AG/DAP | awk '{print $1}')
    nChIP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_AG/ChIP | awk '{print $1}')
    ncommon=$(wc -l $dirout/DAP_${elt}_vs_ChIP_AG/common | awk '{print $1}')
    echo $nDAP $ncommon $nChIP
    python make_venn_figure5.py -DAP $nDAP -ChIP $nChIP -common $ncommon -name_DAP "$(echo -e DAP-Seq\\n${elt})" -name_ChIP "$(echo -e ChIP-Seq\\nAG)" -out $dirout/figures/DAP_${elt}_vs_ChIP_AG.svg -col 8 "${color[${elt}]}"

    ##### plotting reg Venn diagrams #####

    echo "============== now working on regulated genes ======================"

    ##### creating directories #####
    mkdir -p $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg  $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg  $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg  

    # 1 AG reg DAP_${elt}_vs_ChIP_SEP3AG
    grep . $AG_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg/common
    grep . $AG_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/DAP > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg/DAP
    grep . $AG_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/ChIP > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg/ChIP
    nDAP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg/DAP | awk '{print $1}')
    nChIP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg/ChIP | awk '{print $1}')
    ncommon=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_AG_reg/common | awk '{print $1}')
    echo $nDAP $ncommon $nChIP
    python make_venn_figure5.py -DAP $nDAP -ChIP $nChIP -common $ncommon -name_DAP "$(echo -e DAP-Seq\\n${elt})" -name_ChIP "$(echo -e ChIP-Seq\\nSEP3/AG)" -out $dirout/figures/DEG_AG_DAP_${elt}_vs_ChIP_SEP3AG.svg -col 5 "${color[${elt}]}"

    # 2 SEP3 reg DAP_${elt}_vs_ChIP_SEP3AG
    grep . $SEP3_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg/common
    grep . $SEP3_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/DAP > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg/DAP
    grep . $SEP3_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/ChIP > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg/ChIP
    nDAP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg/DAP | awk '{print $1}')
    nChIP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg/ChIP | awk '{print $1}')
    ncommon=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3_reg/common | awk '{print $1}')
    echo $nDAP $ncommon $nChIP
    python make_venn_figure5.py -DAP $nDAP -ChIP $nChIP -common $ncommon -name_DAP "$(echo -e DAP-Seq\\n${elt})" -name_ChIP "$(echo -e ChIP-Seq\\nSEP3/AG)" -out $dirout/figures/DEG_SEP3_DAP_${elt}_vs_ChIP_SEP3AG.svg -col 5 "${color[${elt}]}"

    # 3 SEP3AG reg DAP_${elt}_vs_ChIP_SEP3AG
    grep . $SEP3AG_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/common > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg/common
    grep . $SEP3AG_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/DAP > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg/DAP
    grep . $SEP3AG_reg | xargs -I {} grep {} $dirout/DAP_${elt}_vs_ChIP_SEP3AG/ChIP > $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg/ChIP
    nDAP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg/DAP | awk '{print $1}')
    nChIP=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg/ChIP | awk '{print $1}')
    ncommon=$(wc -l $dirout/DAP_${elt}_vs_ChIP_SEP3AG_SEP3AG_reg/common | awk '{print $1}')
    echo $nDAP $ncommon $nChIP
    python make_venn_figure5.py -DAP $nDAP -ChIP $nChIP -common $ncommon -name_DAP "$(echo -e DAP-Seq\\\n${elt})" -name_ChIP "$(echo -e ChIP-Seq\\nSEP3/AG)" -out $dirout/figures/DEG_SEP3AG_DAP_${elt}_vs_ChIP_SEP3AG.svg -col 5 "${color[${elt}]}"

done


######## make png ###############
mkdir -p $dirout/png

ls $dirout/figures | grep svg | xargs -I {} basename {} .svg | xargs -I {} inkscape -z $dirout/figures/{}.svg -e $dirout/png/{}.png -w 300

exit 0
