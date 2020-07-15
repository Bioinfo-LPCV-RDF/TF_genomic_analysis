#!/bin/bash

meme_matrix=$1
name_matrix=$2

header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE"

start_line=$(grep -n "letter-probability matrix" $meme_matrix | awk -v FS=":" '{print $1}')
matrix=$(tail -n +$start_line $meme_matrix)
line1=$(head -1 <(echo -e "$matrix"))
line1=( $line1 )
nsites=${line1[7]}
pfm=$(tail -n +2 <(echo -e "$matrix") | sed -e 's/^\s\s*//' -e 's/\s\s*/\t/g' -e 's/\s\s*$//g' | awk -v n=$nsites -v OFS="\t" '{print int($1*n)+1,int($2*n)+1,int($3*n)+1,int($4*n)+1}'| sed "1iA\tC\tG\tT")
echo -e "$header\n$pfm"

exit 0

