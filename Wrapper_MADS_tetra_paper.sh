
source /home/304.6-RDF/scripts/DAP_global_analysis_p3.7/compil_functions.sh



main_dir=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7
dir_fastq=$main_dir/Fastq
dir_mapping=$main_dir/Mapping
dir_peakcalling=$main_dir/PeakCalling

############################################################# MAPPING
seed=125
processor=15
fastq_dir="/home/304.3-STRUCTPDEV/shared_data"

fastq1="${fastq_dir}/LIB9-1_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-1_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3rep1 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-3_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-3_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3rep2 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-5_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-5_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3rep3 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB1-5_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB1-5_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3rep4 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-12_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-12_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3AGrep2 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-14_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-14_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3AGrep3 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-8_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-8_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3AGrep1 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-9_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-9_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3delAGrep1 -s $seed -pr $processor

fastq1="${fastq_dir}/LIB9-10_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB9-10_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n SEP3delAGrep2 -s $seed -pr $processor
			
fastq1="${fastq_dir}/LIB1-1_R1_001.fastq.gz"
fastq2="${fastq_dir}/LIB1-1_R2_001.fastq.gz"
main_mapping_Fastq -f1 $fastq1 -f2 $fastq2 -fd $dir_fastq -md $dir_mapping -n control -s $seed -pr $processor

############################################################# PEAKCALLING
Genome_length=120000000
additionnal_args="-q 0.0001 --call-summits"

bam_dir=("$dir_mapping/SEP3rep1" "$dir_mapping/SEP3rep2" "$dir_mapping/SEP3rep3" "$dir_mapping/SEP3rep4")
input=("$dir_mapping/control")
name=SEP3
main_peakcalling -id bam_dir[@] -cd $input -od $dir_peakcalling -nc $name -g $Genome_length -add "$additionnal_args"

bam_dir=("$dir_mapping/SEP3AGrep1" "$dir_mapping/SEP3AGrep2" "$dir_mapping/SEP3AGrep3")
input=("$dir_mapping/control")
name=SEP3AG
main_peakcalling -id bam_dir[@] -cd $input -od $dir_peakcalling -nc $name -g $Genome_length -add "$additionnal_args"

bam_dir=("$dir_mapping/SEP3delAGrep1" "$dir_mapping/SEP3delAGrep2")
input=("$dir_mapping/control")
name=SEP3delAG
main_peakcalling -id bam_dir[@] -cd $input -od $dir_peakcalling -nc $name -g $Genome_length -add "$additionnal_args"

#############################################################  COMPARISONS

initial_comparison -n1 SEP3AG -n2 SEP3delAG -od $main_dir/comparisons -id $dir_peakcalling -r1 2 -r2 0.5
initial_comparison -n1 SEP3AG -n2 SEP3 -od $main_dir/comparisons -id $dir_peakcalling -r1 2 -r2 0.5

#############################################################  Motifs generations (PWM, TFFM)

genome=/home/304.6-RDF/data/tair10.fas
out_motif=${main_dir}/motifs
name=SEP3
peaks_SEP3=$dir_peakcalling/${name}/${name}_narrow.bed
compute_motif -p $peaks_SEP3 -n $name -g $genome -od $out_motif  -ls 600 -nm 1 -mim 16 -mam 16 -s 1

name=SEP3AG
peaks_SEP3AG=$dir_peakcalling/${name}/${name}_narrow.bed
compute_motif -p $peaks_SEP3AG -n $name -g $genome -od $out_motif  -ls 600 -nm 1 -mim 16 -mam 16 -s 1

name=SEP3delAG
peaks_SEP3delAG=$dir_peakcalling/${name}/${name}_narrow.bed
compute_motif -p $peaks_SEP3delAG -n $name -g $genome -od $out_motif  -ls 600 -nm 1 -mim 16 -mam 16 -s 1

name=SEP3AGspe
peaks_SEP3AGspe=$main_dir/comparisons/SEP3AG_SEP3delAG/SEP3AG_spePeaks_2.bed
compute_motif -p $peaks_SEP3AG -n $name -g $genome -od $out_motif  -ls 600 -nm 1 -mim 16 -mam 16 -s 1

#############################################################  Negative set generation

prep_annotation -p 3000  -g $genome -n /home/304.6-RDF/data/tair10 #Folder data must contain tair10.gff
annotation=/home/304.6-RDF/data/tair10.bed

compute_NS -p $peaks_SEP3AG -n SEP3AG -g $genome -od ${main_dir}/NSG -nb 6 -anf $annotation 
compute_NS -p $peaks_SEP3delAG -n SEP3delAG -g $genome -od ${main_dir}/NSG -nb 6 -anf $annotation 
compute_NS -p $peaks_SEP3 -n SEP3 -g $genome -od ${main_dir}/NSG -nb 6 -anf $annotation 
compute_NS -p $peaks_SEP3AGspe -n SEP3AGspe -g $genome -od ${main_dir}/NSG -nb 6 -anf $annotation 

#############################################################  Compute ROCS curves & AUC

peaks=("${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed")
neg_set=("${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed")
matrices=("$out_motif/SEP3AG/SEP3AG.pfm" "$out_motif/SEP3AG/SEP3AG_tffm.xml")
names=("SEP3AG" "SEP3AG")
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3AG
# 
peaks=("${main_dir}/NSG/SEP3delAG_pos.bed" "${main_dir}/NSG/SEP3delAG_pos.bed")
neg_set=("${main_dir}/NSG/SEP3delAG_1_neg.bed" "${main_dir}/NSG/SEP3delAG_1_neg.bed")
matrices=("$out_motif/SEP3delAG/SEP3delAG.pfm" "$out_motif/SEP3delAG/SEP3delAG_tffm.xml")
names=("SEP3delAG" "SEP3delAG")
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3delAG

peaks=("${main_dir}/NSG/SEP3_pos.bed" "${main_dir}/NSG/SEP3_pos.bed")
neg_set=("${main_dir}/NSG/SEP3_1_neg.bed" "${main_dir}/NSG/SEP3_1_neg.bed")
matrices=("$out_motif/SEP3/SEP3.pfm" "$out_motif/SEP3/SEP3_tffm.xml")
names=("SEP3" "SEP3")
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3

peaks=("${main_dir}/NSG/SEP3AGspe_pos.bed" "${main_dir}/NSG/SEP3AGspe_pos.bed")
neg_set=("${main_dir}/NSG/SEP3AGspe_1_neg.bed" "${main_dir}/NSG/SEP3AGspe_1_neg.bed")
matrices=("$out_motif/SEP3AGspe/SEP3AGspe.pfm" "$out_motif/SEP3AGspe/SEP3AGspe_tffm.xml")
names=("SEP3AGspe" "SEP3AGspe")
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3AGspe

#############################################################  Compute ROCS curves & AUC SWAP

# SEP3AG & SEP3delAG matrices on SEP3AG dataset
peaks=("${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed")
neg_set=("${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed")
matrices=("$out_motif/SEP3AG/SEP3AG.pfm" "$out_motif/SEP3AG/SEP3AG_tffm.xml" "$out_motif/SEP3delAG/SEP3delAG.pfm" "$out_motif/SEP3delAG/SEP3delAG_tffm.xml")
names=("SEP3AG" "SEP3AG" "SEP3delAG" "SEP3delAG")
colors=('#40A5C7' '#307C95' '#F9626E' '#BB4A52')
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3delAG_SEP3AG_SEP3AGset -color colors[@]

# SEP3delAG & SEP3AG matrices on SEP3AG dataset
peaks=("${main_dir}/NSG/SEP3delAG_pos.bed" "${main_dir}/NSG/SEP3delAG_pos.bed" "${main_dir}/NSG/SEP3delAG_pos.bed" "${main_dir}/NSG/SEP3delAG_pos.bed")
neg_set=("${main_dir}/NSG/SEP3delAG_1_neg.bed" "${main_dir}/NSG/SEP3delAG_1_neg.bed" "${main_dir}/NSG/SEP3delAG_1_neg.bed" "${main_dir}/NSG/SEP3delAG_1_neg.bed")
matrices=("$out_motif/SEP3delAG/SEP3delAG.pfm" "$out_motif/SEP3delAG/SEP3delAG_tffm.xml" "$out_motif/SEP3AG/SEP3AG.pfm" "$out_motif/SEP3AG/SEP3AG_tffm.xml")
names=("SEP3delAG" "SEP3delAG" "SEP3AG" "SEP3AG")
colors=('#F9626E' '#BB4A52' '#40A5C7' '#307C95')
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3delAG_SEP3AG_SEP3delAGset -color colors[@]

# SEP3AG & SEP3 matrices on SEP3AG dataset
peaks=("${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3AG_pos.bed")
neg_set=("${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3AG_1_neg.bed")
matrices=("$out_motif/SEP3AG/SEP3AG.pfm" "$out_motif/SEP3AG/SEP3AG_tffm.xml" "$out_motif/SEP3/SEP3.pfm" "$out_motif/SEP3/SEP3_tffm.xml")
names=("SEP3AG" "SEP3AG" "SEP3" "SEP3")
colors=('#40A5C7' '#307C95' '#F0875A' '#B46544')
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3_SEP3AG_SEP3AGset -color colors[@]

# SEP3 & SEP3AG matrices on SEP3 dataset
peaks=("${main_dir}/NSG/SEP3_pos.bed" "${main_dir}/NSG/SEP3_pos.bed" "${main_dir}/NSG/SEP3_pos.bed" "${main_dir}/NSG/SEP3_pos.bed")
neg_set=("${main_dir}/NSG/SEP3_1_neg.bed" "${main_dir}/NSG/SEP3_1_neg.bed" "${main_dir}/NSG/SEP3_1_neg.bed" "${main_dir}/NSG/SEP3_1_neg.bed")
matrices=("$out_motif/SEP3/SEP3.pfm" "$out_motif/SEP3/SEP3_tffm.xml" "$out_motif/SEP3AG/SEP3AG.pfm" "$out_motif/SEP3AG/SEP3AG_tffm.xml")
names=("SEP3" "SEP3" "SEP3AG" "SEP3AG")
colors=('#F0875A' '#B46544' '#40A5C7' '#307C95')
compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/SEP3_SEP3AG_SEP3set -color colors[@]

# converting SVG ROCS into PNGs
ls ${main_dir}/ROCS/*.svg | xargs -I {} bash -c "basename {}  .svg" | xargs -I {} --max-procs=20 inkscape -z -e ${main_dir}/ROCS/{}.png -w 1800 -h 1080 ${main_dir}/ROCS/{}.svg

#############################################################  SPACING


peaks=("${main_dir}/NSG/SEP3AG_pos.bed" "${main_dir}/NSG/SEP3delAG_pos.bed" "${main_dir}/NSG/SEP3_pos.bed" "${main_dir}/NSG/SEP3AGspe_pos.bed")
neg_set=("${main_dir}/NSG/SEP3AG_1_neg.bed" "${main_dir}/NSG/SEP3delAG_1_neg.bed" "${main_dir}/NSG/SEP3_1_neg.bed" "${main_dir}/NSG/SEP3AGspe_1_neg.bed")
matrices=("$out_motif/SEP3AG/SEP3AG_tffm.xml" "$out_motif/SEP3delAG/SEP3delAG_tffm.xml" "$out_motif/SEP3/SEP3_tffm.xml" "$out_motif/SEP3AGspe/SEP3AGspe_tffm.xml")
names=("SEP3AG" "SEP3delAG" "SEP3" "SEP3AGspe")
th=("0.4" "0.35" "0.3")
compute_space -p peaks[@] -ns neg_set[@] -nb 6 -m matrices[@] -n names[@] -od ${main_dir}/Spacing -th th[@] -maxy 12 -mins 20 -maxs 80 -ol 9 -or 9

#############################################################  HEATMAP

size_file=/home/304.6-RDF/data/tair10.size # file containing size of each chromosome
./bedGraphToBigWig $dir_peakcalling/SEP3AG/SEP3AG_cov.bdg $size_file $dir_peakcalling/SEP3AG.bw
./bedGraphToBigWig $dir_peakcalling/SEP3delAG/SEP3delAG_cov.bdg $size_file $dir_peakcalling/SEP3delAG.bw

python show_reads.py -b $dir_peakcalling -p ${main_dir}/SEP3AG_SEP3delAG/SEP3AG_SEP3delAG_peaks_reprocessed.bed -r ${main_dir}/SEP3AG_SEP3delAG/ -s 3 -n SEP3AG SEP3delAG


#############################################################  EVALUATION IMPACT OF SPACING
mkdir -p $results/reanalysis/spaced2
results_XPspace=$results/reanalysis/spaced2

table=$results/reanalysis/spaced/SEP3delAG_SEP3AG/table_SEP3delAG_SEP3AG.csv
# 
sed '1d' $table | awk -v OFS="\t" '{print $1,$2,$3}'> $results_XPspace/SEP3delAG_SEP3AG.bed
peaks=$results/reanalysis/spaced2/SEP3delAG_SEP3AG.bed

sed '1d' $table | awk -v OFS="\t" '{print $1,$2,$3,$5,$6}' > $results_XPspace/table_SEP3delAG_SEP3AG.bed
Ntable=$results/reanalysis/spaced2/table_SEP3delAG_SEP3AG.bed

bedtools getfasta -fi /home/304.6-RDF/data/tair10.fas -bed $peaks -fo $results_XPspace/SEP3delAG_SEP3AG.fas
peak_file_fas=$results_XPspace/SEP3delAG_SEP3AG.fas
# 
python get_interdistances.py -tffm $results/reanalysis/motifs/SEP3AG/tffm_first_order.xml  -o $results_XPspace/ -n ${name1}_${name2}_spacing -minInter 20 -maxInter 80  -pos $peak_file_fas  -th 0.001 -neg $peak_file_fas -points True -ol 9 -or 9 -one_panel  -write_inter -no_absolute_panel

th1=0.1
continu1=1
echo -e "score1\tscore2\tcorr" > $results_XPspace/test.txt
while [ $continu1 -eq 1 ];
do
	echo -e "th\t$th1"
	sed 's/[:,/]/\t/g' $results_XPspace/SEP3delAG_SEP3AG_spacing.bed | awk -v th1=$th1 -v th2=$th1 -v OFS="\t" '{save="0";for(i=4;i<=NF;i+=3) {j=i+1;k=i+2;if(($j>=th1 && $k>=th2)||($j>=th2 && $k>=th1)){gsub(/^ER|^IR|^DR/,"",$i);if($i==36||$i==37){if(save==0){save+=1}}else{if($i==46||$i==47){if(save==0){save+=1}}else{if($i==55||$i==56||$i==57||$i==58){if(save==0){save+=1}}else{save+=0}}}}};print $1,$2,$3,save}' | sort -k1,1 -k2,2n | awk -v FS="\t" '{print $4}' > $results_XPspace/Spacing${th1}.inter
	paste $Ntable $results_XPspace/Spacing${th1}.inter | sed "1ichr\tbegin\tend\tSEP3delAG\tSEP3AG\tSpacing" >  $results_XPspace/table_${th1}.csv
	th1=$(calc $th1+0.05)
	continu1=$(awk -vn1="$th1" -vn2="0.9" 'BEGIN{print (n1<=n2)?1:0}')
done
Rscript compute_spacingScorev2.r -f $results_XPspace/table_0.3.csv,$results_XPspace/table_0.35.csv,$results_XPspace/table_0.4.csv,$results_XPspace/table_0.45.csv,$results_XPspace/table_0.5.csv -n SEP3delAG,SEP3AG