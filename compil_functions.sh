# list passed as arg
# list=("a" "b")
# test(){ testVar=("${!1}"); echo "${testVar[@]}";  }
# test list[@]

source compil_usages.sh

###### Needed in path
# fastqc #v0.11.7 mapping_Fastq_bowtie2
# NGmerge #v0.2_dev mapping_Fastq_bowtie2
# bedtools #v2.27.1 peakcalling_MACS2
# mspc #v4.0.0 main_peakcalling
# inkscape #v0.92.2 compute_space
# python 2.7.9  #Python 3 version in progress
# R 3.5.0
PATH_TO_COMPIL=/home/304.3-STRUCTPDEV/MADS/mapping_DAP/pipepline_rep/Scripts_PAPER
# mapping_Fastq_bowtie2
PATH_TO_BOWTIE_INDEX="/home/304.3-STRUCTPDEV/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"
PATH_TO_BOWTIE2=$(dirname $(which bowtie2)) # v2.3.4.1
PATH_TO_SAMTOOLS=$(dirname $(which samtools)) #v1.8 (using htslib v1.8)
PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools))) #v2.27.1
# peakcalling_MACS2
PATH_TO_MACS=$(dirname $(readlink -f $(which macs2))) # v2.1.1.20160309
# initial_comparison
export TMPDIR=/nobackup # use a disk with some space (~20Go for safety)
do_plot_quant=$PATH_TO_COMPIL/bin/Hist_cov_gen.r
merge_peaks=$PATH_TO_COMPIL/bin/merge_all_peaks.py
compute_coverage=$PATH_TO_COMPIL/bin/compute_coverage.py
# Compute motif
meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip # v4.12.0
meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme # v4.12.0
meme2pfm=/home/304.6-RDF/scripts/meme2pfm.sh
prepMEMEforPalTFFM=$PATH_TO_COMPIL/bin/prepMEMEforPalTFFM.py
pfmTOtffm=$PATH_TO_COMPIL/bin/get_tffm.py
# prep_annotation
prepare_gff=$PATH_TO_COMPIL/bin/prepare_gff.py
generate_BedFromGff=$PATH_TO_COMPIL/bin/generate_bed_gff.py
bedops=/home/304.3-STRUCTPDEV/bedops/bin/bedops # v2.4.38
bedmap=/home/304.3-STRUCTPDEV/bedops/bin/bedmap # v2.4.38
# compute_NS
negative_set_script=$PATH_TO_COMPIL/bin/negative_set_generator.py
# compute_ROCS
pocc_pfm=$PATH_TO_COMPIL/bin/compute_POcc.py
tffmscores=$PATH_TO_COMPIL/bin/get_best_score_tffm.py
scores_prog=$PATH_TO_COMPIL/bin/scores.py
plot_ROCS_prog=$PATH_TO_COMPIL/bin/plots_ROCS_multiple.py
# compute_space
spacing_mk=$PATH_TO_COMPIL/bin/get_interdistances.py
tffm_all_scores=/home/304.6-RDF/scripts/get_all_score_tffm.py
plot_spacing_v2=/home/304.6-RDF/scripts/DAP_global_analysis/plot_spacings.r
# heatmap_reads
heatmap_mk=$PATH_TO_COMPIL/bin/show_reads.py
bdg_to_bw=$PATH_TO_COMPIL/bin/bedGraphToBigWig
heatmap_mk=$PATH_TO_COMPIL/bin/compute_spacingScorev2.r

# Function to compute simple math with bash script
# usage : calc $M1/$M2 or calc $M1/$M2*M3
calc(){
awk "BEGIN { print "$*" }"; 
}

# used to create delimited string from array
join_by(){
local IFS="$1"; shift; echo "$*"; 
}

# Function to write to log. Note that if "verbose" is set,
# the log messages will be on the std also
logme(){
    if [ ! -z $verbose ]
       then
           echo "${*}" >> $log 2>&1;
    else
           echo "${*}" >> $log
    fi
}


mapping_Fastq_bowtie2 (){
#mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s <INT> -pr <INT>
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-id)
			local in_dir=$2
			 echo "data directory set to: ${2}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-s)
			local seed=$2
			 echo "data directory set to: ${2}";shift 2;;
		-pr)
			local proc=$2
			 echo "thread(s) number set to: ${2}";shift 2;;
		-h)
			usage mapping_Fastq_bowtie2; exit;;
		--help)
			usage mapping_Fastq_bowtie2; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage mapping_Fastq_bowtie2; exit;;
	esac
done
local Errors=0
if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; local proc=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage mapping_Fastq_bowtie2; exit 1; fi


# PATH_TO_BOWTIE_INDEX="/home/304.3-STRUCTPDEV/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"
# PATH_TO_BOWTIE2=$(dirname $(which bowtie2))
# PATH_TO_SAMTOOLS=$(dirname $(which samtools))
# PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools)))

echo "Processing data set $in_dir";
for fastq in $(find $in_dir -name "*.fastq.gz")
do 
	if [ -e $fastq ]
	then
		local fastq_file=${fastq##*/}
		local fastq_filename=${fastq_file%.*.gz}
		local fileIN=(${fastq_filename//_/ })
		local filename=${fileIN[0]} #get the file name
		local order=0
		for elt in ${fileIN[@]}; 
		do # get the order of a paired-ended FASTQ file
			if [ "R1" = "$elt" ] ; then
				local order=1
				break
			elif [ "R2" = "$elt" ] ; then
				local order=2
				break
			fi
		done
		local paired_file=${fileIN[2]%.*} #get last element which is the name of the paired file
		echo "Step 1: alignment using bowtie2";
		#no paired file
		if [ $order -eq 0 ] #[ "$paired_file" == "$(dirname $fastq)" ]
		then
			local out=$out_dir/$filename #$fastq_filename
			if [ ! -e "$out.sam" ]; 
			then
				$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -U $fastq -S $out.sam -p $proc
			fi
		else
		#figure out the order of the pair-ended files
			local paired_file_order=1
			local fastqIsFirst=false
			if [ "$order" == "$paired_file_order" ]
			then
				local paired_file_order=2
				local fastqIsFirst=true
			fi
		#build up the file name of the pair
			local paired_file=${fastq/"R"$order/"R"$paired_file_order}
		#${paired_file%.*}_${paired_file_order}_$(basename $filename).fastq.gz
			echo "The paired file for $fastq is $paired_file"
			local out=$out_dir/$(basename $filename)_${fileIN[2]%.*} #$(dirname $fastq)/$(basename $filename)_${fileIN[2]%.*}
# 			mkdir -p $out
			local R_NGmerge=$out
			fastqc -o $in_dir $fastq $paired_file 
			# check if the paired file exists in the folder or exit with error
			if [ -e ${paired_file} ] && [ ! -e "$out.sam" ]
			then
				local out=$out_dir/$(basename $filename)_${fileIN[2]%.*} #$(dirname $fastq)/$(basename $filename)_${fileIN[2]%.*}
				if $fastqIsFirst
				then 
					NGmerge -a -1 $fastq -2 $paired_file -q 35 -o $R_NGmerge -n 5
					$PATH_TO_BOWTIE2/bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 ${out}_1.fastq.gz -2 ${out}_2.fastq.gz -S $out.sam --dovetail -p $proc 2> ${out}_log.txt
					#echo "$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -1 $fastq -2 $paired_file -S $fastq.sam"
					echo "Processed file ${out}.sam"
				else
					NGmerge -a -2 $fastq -1 $paired_file -q 35 -o $R_NGmerge -n 5
					$PATH_TO_BOWTIE2/bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 ${out}_1.fastq.gz -2 ${out}_2.fastq.gz -S $out.sam --dovetail -p $proc 2> ${out}_log.txt
					#echo "$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -1 $paired_file -2 $fastq -S $fastq.sam"
					echo "Processed file ${out}.sam" 
				fi
			else
				echo "Paired file for $fastq does not exist for data set $(basename $in_dir)"
				continue
			fi
			#remove the fastq files
			rm $fastq
			rm $paired_file
		fi
		date
		echo "Step 2: filtering SAM";
		$PATH_TO_SAMTOOLS/samtools view -Sh $out.sam | \
		grep -e "^@" -e 'XM:i:[012][^0-9]' | awk '$1~/@/ || $5>30 {print $0}' | grep -v "XS:i:" > $out.filtered.sam    # filter out reads with suboptimal scores
		date 
		echo "Step 3: SAM to BAM conversion"
		$PATH_TO_SAMTOOLS/samtools view -Sh -b $out.filtered.sam \
		> $out.filtered.bam
		$PATH_TO_SAMTOOLS/samtools sort -o $out.filtered.sorted.bam \
		$out.filtered.bam
		date
		echo "Step 4: remove PCR duplicates";
		$PATH_TO_SAMTOOLS/samtools rmdup -s $out.filtered.sorted.bam \
		$out.filtered.sorted.nodup.bam
		date
		echo "Step 5: BAM indexing";
		$PATH_TO_SAMTOOLS/samtools index $out.filtered.sorted.nodup.bam
		date
		echo "Step 6: removing temporary sam"
		rm $out.filtered.sam
		rm $out.sam
		rm ${out}_1.fastq.gz
		rm ${out}_2.fastq.gz
		rm ${out}.filtered.bam
		rm ${out}.filtered.sorted.bam
		date
	fi
done

}


main_mapping_Fastq(){
# main_mapping_Fastq -fd <PATH> -md <PATH> -f1 <FILE> -f2 [FILE] -n <STRING> -s [INT] -pr [INT]
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-fd)
			local Fastq_dir=$2
			 echo "data directory set to: ${2}";shift 2;;
		-md)
			local Mapping_dir=$2
			 echo "Mapping directory set to: ${2}";shift 2;;
		-f1)
			local Fastq1=$2
			 echo "Fastq pair 1 / Fastq non paired end set to: ${2}";shift 2;;
		-f2)
			local Fastq2=$2
			 echo "Fastq pair 2 set to: ${2}";shift 2;;
		-n)
			local Name=$2
			 echo "name of Fastq directory set to: ${2}";shift 2;;
		-s)
			local seed=$2
			 echo "seed for random set to: ${2}";shift 2;;
		-pr)
			local proc=$2
			 echo "thread(s) number set to: ${2}";shift 2;;
		-h)
			usage main_mapping_Fastq; exit;;
		--help)
			usage main_mapping_Fastq; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage main_mapping_Fastq; exit;;
	esac
done
local Errors=0
if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; local proc=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $Fastq2 ]; then echo "-f2 argument not used, assuming non paired-end analysis"; local Fastq2=""; fi
if [ -z $Fastq_dir ]; then echo "ERROR: -fd argument needed"; Errors+=1; fi
if [ -z $Mapping_dir ]; then echo "ERROR: -md argument needed"; Errors+=1; fi
if [ -z $Fastq1 ]; then echo "ERROR: -f1 argument needed"; Errors+=1; fi
if [ -z $Name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage main_mapping_Fastq; exit 1; fi

if [ -d $Mapping_dir/$Name ]; then
  rm -Rf $Mapping_dir/$Name
fi
mkdir -p $Mapping_dir/$Name $Fastq_dir/$Name
if [ $(dirname $Fastq1) != $Fastq_dir/$Name ];then
  cp $Fastq1 $Fastq_dir/$Name
fi
if [ $Fastq2 != "" ]; then
	if [ $(dirname $Fastq2) != $Fastq_dir/$Name ];then
  		cp $Fastq2 $Fastq_dir/$Name
	fi
fi
#mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s <INT> -pr <INT>
mapping_Fastq_bowtie2 -id $Fastq_dir/$Name -od $Mapping_dir/$Name -s $seed -pr $proc

}


peakcalling_MACS2 (){
# peakcalling_MACS2 -id <PATH> -od <PATH> -g <FILE> -add <STRING>
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-id)
			local in_dir=$2
			 echo "data directory set to: ${2}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-g)
			local Genome_length=$2
			 echo "Genome length (mappable) set to: ${2}";shift 2;;
		-add)
			local additionnal_args="$2"
			 echo "additionnal arguments for MACS2: ${2}";shift 2;;
		-h)
			usage peakcalling_MACS2; exit;;
		--help)
			usage peakcalling_MACS2; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage peakcalling_MACS2; exit;;
	esac
done
local Errors=0
if [ -z $Genome_length ]; then echo "-g argument not used, assuming A.thaliana is used: 120000000"; local proc=120000000; fi
if [ -z "$additionnal_args" ]; then echo "-add argument not used, no additional argument passed to MACS2"; local additionnal_args=""; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage peakcalling_MACS2; exit 1; fi

# PATH_TO_MACS=$(dirname $(readlink -f $(which macs2)))
# PATH_TO_SAMTOOLS=$(dirname $(which samtools))
# PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools)))

# START PROCESSING
echo "Processing folder: ${in_dir##*/}"
# add current folder to the output path
local out_dir=${out_dir}/${in_dir##*/}
#check if the output folder is not already there 
if [ ! -d "$out_dir" ]
     then
        #create folder structure
        mkdir -p $out_dir
        mkdir $out_dir/controls
        mkdir $out_dir/replicates     
    
        #split controls and replicates separately (ENCODE datasets have control in the name of the file)     
        for bam in $(find $in_dir -name "*.filtered.sorted.nodup.bam")
           do
             if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" == *"control"* ]]               
                then
                    cp $bam $out_dir/controls
                else
                    cp $bam $out_dir/replicates

		    local libsize=$(samtools view -c $bam)
		    local scale=$(calc 1000000/$libsize)
# 		    bedtools genomecov -bga -scale $scale -ibam $bam -g /home/304.6-RDF/data/tair10.size > $out_dir/${in_dir##*/}_cpm.bdg
		    bedtools genomecov -bga -scale $scale -ibam $bam > $out_dir/${in_dir##*/}_cpm.bdg
            fi              
        done

        local log=$out_dir/${in_dir##*/}.log;

        #check for controls and unzip them or log "no controls"
        local all_ok=true;
        cd $out_dir/controls
        local controls=(*)
        if [ ${#controls[@]} -eq 0 ]
            then
               logme "No control files present"
               local all_ok=false;
            else
#               gunzip *.gz
               local controls=(*)
               #prefix with the absolute path
               local controls=("${controls[@]/#/$out_dir/controls/}")
               local controls=`join_by " " "${controls[@]}"`

               #merge all control files
               $PATH_TO_SAMTOOLS/samtools merge control.bam $controls
               #convert to BED
               $PATH_TO_BEDTOOLS2/bamToBed -i control.bam > control.bed

               #clean up
               rm *.bam
               $PATH_TO_BEDTOOLS2/sortBed -i control.bed > control.bed.sorted
               mv control.bed.sorted control.bed
        fi

        #check for replicates and unzip them or log "no replicates" 
        cd $out_dir/replicates
        local replicates=(*)
        if [ ${#replicates[@]} -eq 0 ]
            then
               logme "No replicates present"
               local all_ok=false;
            else
#               gunzip *.gz;
               local replicates=(*)
               #prefix with the absolute path       
               local replicates=("${replicates[@]/#/$out_dir/replicates/}")
               local replicates=`join_by " " "${replicates[@]}"`
               #merge all replicates files
               $PATH_TO_SAMTOOLS/samtools merge replicate.bam $replicates
               #convert to BED
               $PATH_TO_BEDTOOLS2/bamToBed -i replicate.bam > replicate.bed
               #clean up
               rm *.bam
               $PATH_TO_BEDTOOLS2/sortBed -i replicate.bed > replicate.bed.sorted
               mv replicate.bed.sorted replicate.bed
        fi

        # launch the peak calling by peak caller type
        if [ $all_ok ]
            then
		cd $out_dir
		logme "Started MACS..."
		time $PATH_TO_MACS/macs2 callpeak -t $out_dir/replicates/replicate.bed -c $out_dir/controls/control.bed -B -f BED -n ${in_dir##*/} -g $Genome_length $additionnal_args >> $log 2>&1;
       else
           logme "Not all needed input present."
           exit 1
       fi 
       
       # filter out peaks that are present in Input
       for bdg in $(find $out_dir -name "*_treat_pileup.bdg")
         do
           echo ${in_dir##*/}
           local basename=${bdg%_treat_pileup.bdg}
           local mean_cov=$(awk '{S+=$4}END{print (S/NR)*2}' ${basename}_control_lambda.bdg)
           echo ${basename#${out_dir}/} $mean_cov
           sort -k1,1 -k2,2n ${basename}_peaks.narrowPeak > ${basename}_sorted.narrowPeak
           bedtools intersect -a ${basename}_sorted.narrowPeak -b ${basename}_treat_pileup.bdg -wb | awk -v OFS="\t" '{print $11,$12,$13,$14}' > ${basename}_treatment.bdg
           bedtools intersect -a ${basename}_sorted.narrowPeak -b ${basename}_control_lambda.bdg -wb | awk -v OFS="\t" '{print $11,$12,$13,$14}' > ${basename}_control.bdg
           bedtools unionbedg -i ${basename}_treatment.bdg ${basename}_control.bdg | awk -v OFS="\t" -v Mean=${mean_cov} '$5 > Mean && $2 < $3 && $4 < $5 && $5 >= 1 {print $0}' > ${basename}_forbidden_regions.bed
           bedtools intersect -v -a ${basename}_peaks.narrowPeak -b ${basename}_forbidden_regions.bed -wa > ${basename}_filtered.narrowPeak
           bedtools intersect -a ${basename}_peaks.narrowPeak -b ${basename}_forbidden_regions.bed -wa  > ${basename}_removed.narrowPeak
         done
       
#       clean up
       rm -R $out_dir/controls;
       rm -R $out_dir/replicates;

       awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' ${basename}_filtered.narrowPeak > ${basename}_peaks.bed
else
    echo "Folder $in_dir already exists in results. Skipping...";
fi

}


main_peakcalling (){
# main_peakcalling -id -cd -od -nc -g -add
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-id)
			local in_dir=("${!2}")
			 echo "data directory set to: ${2}";shift 2;;
		-cd)
			local input_dir=$2
			 echo "control directory set to: ${2}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-nc)
			local name_cons=$2
			 echo "name for consensus directory set to: ${2}";shift 2;;
		-g)
			local Genome_length=$2
			 echo "genome length (mappable) set to: ${2}";shift 2;;
		-add)
			local additionnal_args="$2"
			 echo "additionnal argument set to MACS2: ${2}";shift 2;;
		-top)
			local top=$2
			echo "maximum number of peaks set to: ${2}";shift 2;;
		-h)
			usage main_peakcalling; exit;;
		--help)
			usage main_peakcalling; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage main_peakcalling; exit;;
	esac
done
local Errors=0
if [ -z $Genome_length ]; then echo "-g argument not used, assuming A.thaliana is used: 120000000"; local proc=120000000; fi
if [ -z "$additionnal_args" ]; then echo "-add argument not used, no additional argument passed to MACS2"; local additionnal_args=""; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $input_dir ]; then echo "ERROR: -cd argument needed"; Errors+=1; fi
if [ -z $name_cons ]; then echo "ERROR: -nc argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage main_peakcalling; exit 1; fi
local list_peaks_mspc=()
local list_bdg=()

for dir in ${in_dir[@]};
do
	bam_ech=$(find $dir -name "*.filtered.sorted.nodup.bam")
  for bam in $(find $input_dir -name "*.filtered.sorted.nodup.bam")
  do
    cp $bam $dir/control.filtered.sorted.nodup.bam
  done
  local name=${dir##*/}
  list_peaks_mspc+=("-i $out_dir/$name/${name}_peaks.bed")
  list_bdg+=("$out_dir/$name/${name}_cpm.bdg")
  peakcalling_MACS2 -id $dir -od $out_dir -g $Genome_length -add "$additionnal_args"

for elt in ${bam_ech[@]};
do
	if [ $elt != ${dir}/control.filtered.sorted.nodup.bam ]; then
		total=$(samtools view -c $elt)
		inpeak=$(bedtools sort -i $out_dir/$name/${name}_peaks.bed | bedtools merge -i stdin | bedtools intersect -u -a $elt -b stdin -ubam | samtools view -c)
		FRIP=$(calc $inpeak/$total*100)
		echo "$name $FRIP% $inpeak / $total" > $out_dir/$name/${name}_FreqReadInPeak.txt
	fi
done

done
local list_peaks_mspc=`join_by " " "${list_peaks_mspc[@]}"`

if [[ $(uname -a) =~ "el7" ]]; then
  mspc ${list_peaks_mspc} -r Tec -w 1e-4 -s 1e-8 -c 100% -o $out_dir/$name_cons
else
  echo "Error SL7 node is needed for mspc" && exit 1
fi
if [ -z $top ]; then  local top=$(wc -l $out_dir/$name_cons/ConsensusPeaks.bed) ; fi
sed '1d' $out_dir/$name_cons/ConsensusPeaks.bed | sort -k5,5nr | awk -v OFS="\t" '{print $1,$2,$3}' > $out_dir/$name_cons/${name_cons}_comp.bed
head -${top} $out_dir/$name_cons/${name_cons}_comp.bed > $out_dir/$name_cons/$name_cons.bed

mkdir -p $out_dir/$name_cons/temp
local i=0
local files=()
for dir in ${in_dir[@]};
do
  local name=${dir##*/}
  if [ $i -eq 0 ];
  then
    bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4}END{print save}' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed 
  else
    bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4}END{print save}' | awk -v OFS=" " '{$1=$2=$3="";print $0}' | sed 's/   //' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed
  fi
  files+=("$out_dir/$name_cons/temp/tmp_peaks_$i.bed")
  local i=$(($i+1))
done

paste ${files[@]} > $out_dir/$name_cons/temp/${name_cons}_narrow.bed
awk -v OFS="\t" 'function abs(v) {return v < 0 ? -v : v} {nbmax=0; for(i=4;i<=NF;i++) { if(nbmax==0){max[i-3]=$i; moy[$i]=$i;nb[$i]=1;nbmax=1} if(nbmax!=0){ ok=1;for(lmax in max){ localmax=max[lmax]; if(abs(localmax-$i)<=200){ moy[localmax]+=$i; nb[localmax]+=1; ok=0; break }} if(ok==1){ max[i-3]=$i; moy[$i]=$i;nb[$i]=1}  }}; for(kmax in max){localmax=max[kmax];print($1, int(moy[localmax]/nb[localmax])-200, int(moy[localmax]/nb[localmax])+200); delete max[kmax];delete nb[localmax]; delete moy[localmax]}}' $out_dir/$name_cons/temp/${name_cons}_narrow.bed | sed 's/\t\+/\t/g;s/^\t//' > $out_dir/$name_cons/${name_cons}_narrow.bed
local list_bdg=`join_by " " "${list_bdg[@]}"`
bedtools unionbedg -i ${list_bdg} > $out_dir/$name_cons/${name_cons}_cov_sep.bdg
awk -v OFS="\t" '{moy=0;nb=0;for(i=4;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/$name_cons/${name_cons}_cov_sep.bdg > $out_dir/$name_cons/${name_cons}_cov.bdg
echo "here2"
local i=0
local files=()
for dir in ${in_dir[@]};
do
  local name=${dir##*/}
  if [ $i -eq 0 ];
  then
	bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cpm.bdg -wao | awk -v OFS="\t" '{print $1":"$2"-"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 > $out_dir/$name_cons/temp/tmp_max_${i}.txt
  else
    bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cpm.bdg -wao | awk -v OFS="\t" '{print $1":"$2"-"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 | awk '{print $2}' > $out_dir/$name_cons/temp/tmp_max_${i}.txt
  fi 
  files+=("$out_dir/$name_cons/temp/tmp_max_${i}.txt")
  local i=$(($i+1))
done
paste ${files[@]} | sed "s/[-:]/\t/g" > $out_dir/$name_cons/${name_cons}_max.bed
}


initial_comparison(){
# initial_comparison -n1 -n2 -od -id -r1 -r2
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n1)
			local name1=$2
			echo "name of directory for dataset 1 is: ${2}";shift 2;;
		-n2)
			local name2=$2
			echo "name of directory for dataset 2 is: ${2}";shift 2;;
		-od)
			local result=$2
			echo "output directory set to: ${2}";shift 2;;
		-id)
			local data=$2
			echo "general data directory set to: ${2}";shift 2;;
		-r1)
			local ratio1=$2
			echo "ratio of coverage for specific peaks of dataset 1 is: ${2}";shift 2;;
		-r2)
			local ratio2=$2
			echo "ratio of coverage for specific peaks of dataset 2 is: ${2}";shift 2;;
		-f)
			local filterCov=$2
			echo "coverage must be at least >${2} in both sample to be considered as peaks";shift 2;;
		-he)
			local filterHeight=$2
			echo "height must be at least >${2} in both sample to be considered as peaks";shift 2;;
		-h)
			usage initial_comparison; exit;;
		--help)
			usage initial_comparison; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage initial_comparison; exit;;
	esac
done
local Errors=0
if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
if [ -z $result ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $data ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $filterCov ]; then echo "no filter on coverage applied"; local filterCov=0; fi
if [ -z $filterHeight ]; then echo "no filter on peak height applied"; local filterHeight=0; fi
if [ $Errors -gt 0 ]; then usage initial_comparison; exit 1; fi
if [ -z $ratio1 ]; then 
	if [ -z $ratio2 ]; then 
		echo "no ratios specified, using 2 fold for both"
		local ratio1=2
		local ratio2=0.5
	else
		echo "WARNING: only ratio for $name2 specified, assuming equivalent fold ratio for $name1"
		local ratio1=$(calc 1/${ratio2})
	fi
fi
if [ -z $ratio2 ]; then 
	echo "WARNING: only ratio for $name1 specified, assuming equivalent fold ratio for $name2"
	local ratio2=$(calc 1/${ratio1})
fi

# 	export TMPDIR=/nobackup
# 	local do_plot_quant=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/Hist_cov_gen.r
# 	local merge_peaks=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/merge_all_peaks.py
# 	local compute_coverage=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/compute_coverage.py

	local out_dir=$result/${name1}_${name2}
	mkdir -p $out_dir $result/plots
	local cov1=$data/$name1/${name1}_cov.bdg
	local cov2=$data/$name2/${name2}_cov.bdg

	###### few lines added for selecting peaks on their height:
	if [ $filterHeight -eq 0 ];then
			local peaks1=$data/$name1/${name1}_narrow.bed
    		local peaks2=$data/$name2/${name2}_narrow.bed
	else
# 			awk -v filt=$filterHeight -v OFS="\t" '$4>filt && $5>filt && $6>filt' $data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
			local peaks1=$data/$name1/${name1}_narrow.heightFiltered.bed
# 			awk -v filt=$filterHeight -v OFS="\t" '$4>filt && $5>filt && $6>filt' $data/$name2/${name2}_max.bed > 	$data/$name2/${name2}_narrow.heightFiltered.bed
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/$name2/${name2}_max.bed > $data/$name2/${name2}_narrow.heightFiltered.bed
			local peaks2=$data/$name2/${name2}_narrow.heightFiltered.bed
        fi

    awk -v OFS="\t" -v name=$name1 '{print $1,$2,$3,name}' $peaks1  > $out_dir/${name1}_peaks.bed
    awk -v OFS="\t" -v name=$name2 '{print $1,$2,$3,name}' $peaks2  > $out_dir/${name2}_peaks.bed
    
    cat $out_dir/${name1}_peaks.bed $out_dir/${name2}_peaks.bed | sort -k1,1 -k2,2n > $out_dir/${name1}_${name2}_peaks.bed
    python $merge_peaks -f1 $name1 -f2 $name2 -o $result
    echo "merged"
    local peak_file=$out_dir/${name1}_${name2}_peaks_processed.bed
    cat $out_dir/${name1}_${name2}_peaks_merged.bed $out_dir/${name1}_${name2}_peaks_uniques.bed | sort -k1,1 -k2,2n > $peak_file
    echo "sorted"
    bedtools intersect -a $peak_file -b $cov1 -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/$name1.inter &
    bedtools intersect -a $peak_file -b $cov2 -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/$name2.inter &
    wait
    echo "intersected"
    python $compute_coverage -i $out_dir/$name1.inter -m 
    python $compute_coverage -i $out_dir/$name2.inter -m &
    wait 
    echo "cov computed"
    awk  '{print $5/($3-$2)}' $out_dir/$name1.inter.cov  > $out_dir/$name1.inter.cov.normed # adjust cov by len peak
    awk  '{print $5/($3-$2)}' $out_dir/$name2.inter.cov  > $out_dir/$name2.inter.cov.normed
    echo "norm adjusted"
	if [ $filterCov -eq 0 ];then
		paste $peak_file $out_dir/$name1.inter.cov.normed $out_dir/$name2.inter.cov.normed  |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" >  $out_dir/table_${name1}_${name2}.csv
	else
		paste $peak_file $out_dir/$name1.inter.cov.normed $out_dir/$name2.inter.cov.normed | awk -v filter=$filterCov -v OFS="\t" '$5>=filter || $6>=filter {print $0}' |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" >  $out_dir/table_${name1}_${name2}.csv
	fi
    Rscript $do_plot_quant $result $name1 $name2 $out_dir/table_${name1}_${name2}.csv $ratio1 $ratio2
}


compute_motif(){

local pal=false
# compute_motif -p -n -g -od -ls -nm -mim -mam -pal -s
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			echo "peak file set to: ${2}";shift 2;;
		-n)
			local name=$2
			echo "name for sub-directory set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "genome to use set to: ${2}";shift 2;;
		-od)
			local results_motif=$2
			echo "result directory set to: ${2}";shift 2;;
		-ls)
			local learning_size=$2
			echo "size for learning set: ${2}";shift 2;;
		-nm)
			local nb_motifs=$2
			echo "number of PWM motif generated: ${2}";shift 2;;
		-mt)
			local motifstffm=$2
			echo "PWM motif used for TFFM generation set to: ${2}";shift 2;;
		-mim)
			local min_motif=$2
			echo "minimum length for motif set to: ${2}";shift 2;;
		-mam)
			local max_motif=$2
			echo "maximum length for motif set to: ${2}";shift 2;;
		-pal)
			local pal=true
			echo "palindromic mode activated";shift 1;;
		-s)
			local seed=$2
			echo "seed for random set to: ${2}";shift 2;;
		-h)
			usage compute_motif; exit;;
		--help)
			usage compute_motif; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_motif; exit;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $results_motif ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $learning_size ]; then echo "-ls argument not used, using 600 peaks"; local  learning_size=600; fi
if [ -z $motifstffm ]; then echo "-mt argument not used, first motif used for tffm creation"; local motifstffm=0; fi
if [ -z $min_motif ]; then echo "-mim argument not used, min size of motif is 6"; local min_motif=6; fi
if [ -z $max_motif ]; then echo "-mam argument not used, max size of motif is 15"; local min_motif=15; fi
if [ -z $nb_motifs ]; then echo "-nm argument not used, searching for 1 motif"; local nb_motifs=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if ! $pal ; then echo "palindromic mode not activated"; fi
if [ $Errors -gt 0 ]; then usage prep_annotation compute_motif; exit 1; fi

if [ $max_motif -lt $min_motif ]; then
	echo "WARNING: things got mixed up, max size of motif can't be lower than min size."
	echo -e "\tUsing min size as max size "
	local max_motif=$(calc $max_motif+$min_motif) # max motif = total
	local min_motif=$(calc $max_motif-$min_motif) # new min_motif = old max_motif (total - old min = old max)
	local max_motif=$(calc $max_motif-$min_motif) # new max_motif = old min_motif (total - new min = old min)
fi

# prog used here
# meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip
# meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme
# meme2pfm=/home/304.6-RDF/scripts/meme2pfm.sh
# # meme2pfm=/home/304.6-RDF/scripts/DAP_global_analysis/meme2pfm_modifJM.sh
# prepMEMEforPalTFFM=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/prepMEMEforPalTFFM.py
# pfmTOtffm=/home/304.6-RDF/LFY/scripts/get_tffm.py

mkdir -p $results_motif/$name/sets $results_motif/$name/meme $results_motif/$name/tffm
    
    head -$learning_size $peaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
	sed "1,${learning_size}d" $peaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
    bedtools getfasta -fi $genome -fo $results_motif/$name/sets/${name}.fas -bed $results_motif/$name/sets/${name}.bed
    
    if $pal ; then
    echo "pal"
    $meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-pal -meme-maxsize $(calc $learning_size*1000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed 
    else
    $meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-maxsize $(calc $learning_size*1000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed
#     $meme_prog -oc $results_motif/$name/meme -meme-searchsize $learning_size -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs 1 -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed
    fi
    $meme2meme $results_motif/$name/meme/meme_out/meme.txt > $results_motif/$name/meme/meme_out/meme_mini.txt

local path_to_meme_mini=$results_motif/$name/meme/meme_out/meme_mini.txt

bash $meme2pfm $path_to_meme_mini $name > $results_motif/$name/${name}.pfm

if $pal ; then
  python $prepMEMEforPalTFFM -m $results_motif/$name/meme/meme_out/meme.txt -f $results_motif/$name/sets/${name}_tffm_learningset.fas
  local learning_set_tffm=$results_motif/$name/sets/${name}_tffm_learningset.fas
else
  local learning_set_tffm=$results_motif/$name/sets/${name}.fas
fi

python $pfmTOtffm -r $results_motif/$name/tffm/ -f $learning_set_tffm -m $results_motif/$name/meme/meme_out/meme.txt -p 0

cp $results_motif/$name/meme/meme_out/logo1.png $results_motif/$name/${name}_logo.png
cp $results_motif/$name/tffm/tffm_first_order.xml $results_motif/$name/${name}_tffm.xml
}


prep_annotation(){

while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local bed_genome=$2
			echo "gff file of the genome set to: ${2}";shift 2;;
        -p)
			local promoterLength=$2
			echo "osize for the promoters set to: ${2}";shift 2;;
        -g)
			local Genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
		-h)
			usage prep_annotation; exit;;
		--help)
			usage prep_annotation; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage prep_annotation; exit;;
	esac
done

local Errors=0
if [ -z $bed_genome ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $Genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
if [ -z $promoterLength ]; then echo "-p argument not used, using 1000bp as promoter length";local promoterLength=1000 ;fi
if [ $Errors -gt 0 ]; then usage prep_annotation; exit 1; fi

# prepare_gff=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/prepare_gff.py
# generate_BedFromGff=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/generate_bed_gff.py
# bedops=/home/304.3-STRUCTPDEV/bedops/bin/bedops
# bedmap=/home/304.3-STRUCTPDEV/bedops/bin/bedmap

python $prepare_gff -gi ${bed_genome}.gff3 -go ${bed_genome}_prep.gff3  -f ${Genome}
python $generate_BedFromGff -g ${bed_genome}_prep.gff3 -p $promoterLength -o ${bed_genome}_prep.bed

bedtools sort -i ${bed_genome}_prep.bed | $bedops --partition - | $bedmap --echo --echo-map-id --delim '\t' - ${bed_genome}_prep.bed | awk -v OFS="\t" '$4==""{print $1,$2,$3,"intergenic";next} {print $0}' > ${bed_genome}.bed


}


compute_NS(){
# compute_NS -p -n -g -od -anf -nb -bs -gc -lt -s
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			echo "peaks file set to: ${2}";shift 2;;
		-n)
			local name=$2
			echo "prefix for results set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-anf)
			local annotation_file=$2
			echo "Annotation file set to: ${2}";shift 2;;
		-nb)
			local number_of_NS=$2
			echo "number of negative sets to create set to: ${2}";shift 2;;
		-ws)
			local window_size=$2
			echo "window size set to: ${2}";shift 2;;
		-gc)
			local deltaGC=$2
			echo "delta GC allowed set to: ${2}";shift 2;;
		-lt)
			local limit_type=$2
			echo "number of region by origin set to: ${2}";shift 2;;
		-s)
			local seed=$2
			echo "seed for random set to: ${2}";shift 2;;
		-h)
			usage compute_NS; exit;;
		--help)
			usage compute_NS; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_NS; exit;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $annotation_file ]; then echo "ERROR: -anf argument needed"; Errors+=1; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, searching for 1 negative set"; local number_of_NS=1; fi
if [ -z $window_size ]; then echo "-ws argument not used, using 250 bp as window size"; local window_size=250 ; fi
if [ -z $deltaGC ]; then echo "-gc argument not used, using 0.03 as maximum of GC% divergence"; local deltaGC=0.03 ; fi
if [ -z $limit_type ]; then echo "-lt argument not used, using 1000 as minimum group size for origin (use -h for more information)"; local limit_type=1000 ; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ $Errors -gt 0 ]; then usage compute_NS; exit 1; fi

# negative_set_script=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/negative_set_generator.py

mkdir -p ${results}
# awk -v OFS='\t' ' ($3-$2)>1200 {print $1, ($2+($3-$2)-600), ($2+($3-$2)+600);next} {print $1, $2, $3}' $peak > ${results}/${name}_fil.bed
python $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size

}


compute_ROCS(){

# compute_ROCS -p -ns -m -n -g -od -pc
local pocc=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "fasta of the genome set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-pc)
			local pocc=true
			echo "Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
		-color)
			local colors=("${!2}")
			echo "colors set to: ${2}";shift 2;;
		-h)
			usage compute_ROCS; exit;;
		--help)
			usage compute_ROCS; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_ROCS; exit;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if [ -z $colors ]; then echo "-color argument not used, using default"; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544'); fi

if [ $Errors -gt 0 ]; then usage compute_ROCS; exit 1; fi

# pocc_pfm=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/compute_POcc.py
# tffmscores=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/get_best_score_tffm.py
# scores_prog=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/scores.py
# plot_ROCS_prog=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/plots_ROCS_multiple.py

local scores=()
local endnames=()
i=0
mkdir -p ${results}/scores
# for peak in ${peaks[@]}
for ((i=0;i<${#peaks[@]};i++))
do
	local negative_set=${negative_sets[i]}
	local name=${names[i]} # TODO replace \s in name by _ 
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
    echo -e "\tpeaks $peak"
    echo -e "\tNegset $negative_set"
    echo -e "\tname $name"
    echo -e "\tmatrice $matrice"
	if [[ $peak != *".fa"* ]]; then
		bedtools getfasta -fi $genome -bed $peak -fo $results/pos_set.fa
		local peak=$results/pos_set.fa
		bedtools getfasta -fi $genome -bed $negative_set -fo $results/neg_set.fa
		local negative_set=$results/neg_set.fa
	fi

	if [[ $matrice == *".xml"* ]]; then # TFFM scores computation
		python $tffmscores -o ${results}/scores/tffm_scores_pos.tsv -pos $peak -t ${matrice}
		python $tffmscores -o ${results}/scores/tffm_scores_neg.tsv -pos $negative_set -t ${matrice}
		paste <(awk '{print $8}' ${results}/scores/tffm_scores_pos.tsv) <(awk '{print $8}' ${results}/scores/tffm_scores_neg.tsv) >"${results}/scores/TFFM_${name}_scores.tsv"
		scores+=("${results}/scores/TFFM_${name}_scores.tsv")
		endnames+=("${name}_TFFM")
	fi
	if [[ $matrice == *".txt"* ]]; then # K-mer scores computation
		echo "WIP"
	fi
	if [[ $matrice == *".pfm"* ]]; then  # PWM scores
		python $scores_prog -m ${matrice} -f $peak -o ${results}/scores/
		python $scores_prog -m ${matrice} -f $negative_set -o ${results}/scores/
		paste <(sort -u -k1,1 -k4,4nr ${results}/scores/$(basename $peak).scores | sort -u -k1,1 | awk '{print $4}'  | sort -nr )  <(sort -u -k1,1 -k4,4nr ${results}/scores/$(basename $negative_set).scores | sort -u -k1,1 | awk '{print $4}' | sort -nr ) | awk '{print $0}' >"${results}/scores/tab_${name}.tsv"
		scores+=("${results}/scores/tab_${name}.tsv")
		endnames+=("${name}_PWM")
		if $pocc ; then # PWM Pocc
			python $pocc_pfm -s ${results}/scores/$(basename $peak).scores -o ${results}/scores/Pocc_pos.pocc
			python $pocc_pfm -s ${results}/scores/$(basename $negative_set).scores -o ${results}/scores/Pocc_neg.pocc
			paste <( sort -nr ${results}/scores/Pocc_pos.pocc )  <( sort -nr ${results}/scores/Pocc_neg.pocc ) > "${results}/scores/${name}_Pocc_scores.tsv"
			scores+=("${results}/scores/${name}_Pocc_scores.tsv")
			endnames+=("${name}_Pocc")
		fi
	fi
	
done

local scores=`join_by " " "${scores[@]}"`
local endnames=`join_by " " "${endnames[@]}"`
local list_colors=`join_by " " "${colors[@]}"`
python $plot_ROCS_prog -s $scores -n $endnames -o ${results} -of ROC.svg -c ${list_colors}
}

compute_space(){

# compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -max <MAXY>
local pocc=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${2}";shift 2;;
        -th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${2}";shift 2;;
        -maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "number of negative files to use set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-pc)
			local pocc=true
			echo "Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
		-h)
			usage compute_space; exit;;
		--help)
			usage compute_space; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_space; exit;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $thresholds ]; then echo "ERROR: -th argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; local maxy=0; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; local maxSpace=50; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; local minSpace=0; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; local offset_left=0; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; local offset_right=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, 1 negative set used"; local number_of_NS=1; fi
# if [ -z $colors ]; then echo "-od argument not used, "; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544'); fi TODO implement

if [ $Errors -gt 0 ]; then usage compute_space; exit 1; fi

mkdir -p $results
# spacing_mk=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/get_interdistances.py

for ((i=0;i<${#peaks[@]};i++))
do
	local negative_set=${negative_sets[i]}
	local name=${names[i]} # TODO replace \s in name by _ 
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	
	local th=`join_by " " "${thresholds[@]}"`
	neg_set="NA"
	if [[ $peak != *".fa"* ]]; then
        bedtools getfasta -fi $genome -bed $peak -fo $results/pos_set.fa
        local peak=$results/pos_set.fa
    fi
	for ((j=1;j<$number_of_NS;j++))
	do
        if [[ ${negative_set} != *".fa"* ]];then
            bedtools getfasta -fi $genome -bed ${negative_set/_1_neg.bed/_${j}_neg.bed} -fo $results/neg_set_${j}.fa
            if [[ ${neg_set} == "NA" ]]; then
                neg_set="$results/neg_set_${j}.fa"
            else 
                neg_set+=" $results/neg_set_${j}.fa"
            fi
        else
			if [[ ${neg_set} == "NA" ]]; then
				neg_set="  ${negative_set/_1_neg.fa/_${j}_neg.fa}"
			else
            	neg_set+="  ${negative_set/_1_neg.fa/_${j}_neg.fa}"
			fi
        fi
	done
	echo $matrice
	if [[ $matrice == *".pfm"* ]]; then
        python $spacing_mk -mat $matrice -o $results -n ${4%.*}_pwm.svg -minInter $minSpace -maxInter $maxSpace  -pos $results/sets/interdist_peaks_pos.fas -th $th -neg $neg_set -points True -no_absolute_panel -ol $offset_left -or $offset_right --write_inter -one_panel -maxy $maxy
	fi
	if [[ $matrice == *".xml"* ]]; then
        python $spacing_mk -tffm $matrice -o $results -n ${name}_tffm.svg -minInter $minSpace -maxInter $maxSpace  -pos $peak -th $th -neg $neg_set -points True -no_absolute_panel -ol $offset_left -or $offset_right --write_inter -one_panel -maxy $maxy
    fi
	inkscape -z -e $results/${name}_tffm.png -w 1800 -h 1080 $results/${name}_tffm.svg
done
}


compute_space_v2(){

# compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -maxy <MAXY> -maxs -mins -ol -or -g -wi -nap -op -co
local write_inter=false
local no_absolute_panel=false
local one_panel=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${2}";shift 2;;
		-th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${2}";shift 2;;
		-maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "number of negative files to use set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-wi)
			local write_inter=true
			echo "All spacings will be reported in file";shift 1;;
		-nap)
			local no_absolute_panel=true
			echo "absolute enrichment panel won't be plotted";shift 1;;
		-op)
			local one_panel=true
			echo "relative enrichment panels will be merged";shift 1;;
		-co)
			local colors=("${!2}")
			echo "Colors to use: ${2}"; shift 2;;
		-h)
			usage compute_space_v2; exit;;
		--help)
			usage compute_space_v2; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_space; exit;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $names ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $thresholds ]; then echo "ERROR: -th argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; local maxy=0; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; local maxSpace=50; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; local minSpace=0; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; local offset_left=0; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; local offset_right=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, 1 negative set used"; local number_of_NS=1; fi
if [ -z $colors ]; then echo "-co argument not used, "; local colors=('#004899' '#004088' '#003876' '#0050aa' '#1961b2' '#3272bb' '#4c84c3' '#6696cc' '#7fa7d4' '#99b9dd'); fi

if [ $Errors -gt 0 ]; then usage compute_space_v2; exit 1; fi


local additional_python=""
local additional_R=()
if [ $write_inter = true ];then
	additional_python+=("-wi")
fi
if [ $no_absolute_panel = true ];then
	additional_python+=("-nap")
	additional_R+=("--no_absolute_panel")
fi
if [ $one_panel = true ];then
	additional_R+=("--one_panel")
fi
local list_additional_R=`join_by " " "${additional_R[@]}"`
local list_additional_python=`join_by " " "${additional_python[@]}"`
local list_colors=`join_by "," "${colors[@]}"`

for ((i=0;i<${#peaks[@]};i++))
do
	
	local negative_set=${negative_sets[i]}
	local name=${names[i]// /_}
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	local Neg_is_fasta="false"
	local length_mat=0
	mkdir -p $results/${name}/scores
	
	local th=`join_by " " "${thresholds[@]}"`
	neg_set="NA"
	if [[ $peak == *"score"* ]];then
		local pos_file=$peak # no preparation needed
	elif [[ $peak != *".fa"* ]]; then
		if [ ! -f $results/${name}/${name}_pos_set.fa ];then
			bedtools getfasta -fi $genome -bed $peak -fo $results/${name}/${name}_pos_set.fa # need to tranform in fasta
		fi
		local peak=$results/${name}/${name}_pos_set.fa
	fi
	if [[ $peak == *".fa"* ]]; then # need to compute scores
		if [[ $matrice == *".xml"* ]]; then
			if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_pos.tsv ];then
				python $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $peak -t ${matrice}
			fi
			local pos_file=${results}/${name}/scores/${name}_tffm_scores_pos.tsv
		elif [[ $matrice == *".pfm"* ]]; then
			local length_mat=$(awk -v OFS="[\t ]" '{if(NR==1){if($5=="SIMPLE"){typM="simple"} else {typM="dependency"};count=-1;next};if(typM=="simple"){count++};if(typM=="dependency"){if($1=="DEPENDENCY"){count-=2;exit};count++}}END{print count}')
			if [ ! -f ${results}/${name}/scores/$(basename $peak).scores ];then
				python $scores_prog -m ${matrice} -f $peak -o ${results}/${name}/scores/
			fi
			local pos_file=${results}/${name}/scores/$(basename $peak).scores
		fi
	fi

	for ((j=1;j<=$number_of_NS;j++))
	do
		local tmp_fasta="NA"
		if [[ ${negative_set} == *"score"* ]];then
			local tmp_neg_set="${negative_set/_1_neg.fa.scores/_${j}_neg.fa.scores}"
		elif [[ ${negative_set} == *".bed"* ]];then
			if [ ! -f $results/${name}/${name}_set_${j}_neg.fa ];then
				bedtools getfasta -fi $genome -bed ${negative_set/_1_neg.bed/_${j}_neg.bed} -fo $results/${name}/${name}_set_${j}_neg.fa
			fi
			local tmp_fasta=$results/${name}/${name}_set_${j}_neg.fa
			local Neg_is_fasta="true"
		elif [[ ${negative_set} == *".fa"* ]];then
			local Neg_is_fasta="true"
			local tmp_fasta=${negative_set/_1_neg.fa/_${j}_neg.fa}
		fi
		if [ $Neg_is_fasta == "true" ]; then
			if [[ $matrice == *".xml"* ]]; then
				if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv ];then
					python $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv -pos $tmp_fasta -t ${matrice}
				fi
				local tmp_neg_set="${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv"
			elif [[ $matrice == *".pfm"* ]]; then
				if [ ! -f ${results}/${name}/scores/$(basename $tmp_fasta).scores ];then
					python $scores_prog -m ${matrice} -f $tmp_fasta -o ${results}/scores/
				fi
				local tmp_neg_set="${results}/${name}/scores/$(basename $tmp_fasta).scores"
			fi
		fi
		if [[ ${neg_set} == "NA" ]]; then
			local neg_set="$tmp_neg_set"
		else
			local neg_set+=" $tmp_neg_set" 
		fi
	done 
 	python $spacing_mk -o $results/${name}/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th $th -neg $neg_set -ol $offset_left -or $offset_right -lm $length_mat $list_additional_python
	Rscript $plot_spacing_v2 -e ${results}/${name}/${name}_enrichment.tsv -d ${results}/${name}/${name} -r ${results}/${name}/${name}_rate.tsv -m $minSpace -c $list_colors -y $maxy $list_additional_R
done


}


heatmap_reads(){

# heatmap_reads
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-b)
			local bedgraphs=("${!2}")
			echo "bedgraphs files set to: ${2}";shift 2;;
        -n)
			local names=("${!2}")
			echo "names set to: ${2}";shift 2;;
		-od)
			local out_dir=$2
			echo "output directory set to: ${2}";shift 2;;
		-or)
			local order=$2
			echo "order set to: ${2}";shift 2;;
		-s)
			local size_file=$2
			echo "file of chromosome sizes set to: ${2}";shift 2;;
		-cp)
			local comp_dir=$2
			echo "comparison directory set to: ${2}";shift 2;;
		-ws)
			local window_size=`calc $2/2`
			echo "window size around center of peaks set to: ${2}";shift 2;;
		-h)
			usage heatmap_reads; exit;;
		--help)
			usage heatmap_reads; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage heatmap_reads; exit;;
	esac
done

local Errors=0
if [ -z $bedgraphs ]; then echo "ERROR: -b argument needed";Errors+=1;fi
if [ -z $names ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ -z $size_file ]; then echo "ERROR: -s argument needed";Errors+=1;fi
if [ -z $comp_dir ]; then echo "ERROR: -cp argument needed";Errors+=1;fi

if [ -z $window_size ]; then echo "-ws argument not used, using 1000bp";local window_size=1000 ;fi
if [ -z $order ]; then echo "-or argument not used, ordering by CFR";local order=3 ;fi
if [ $Errors -gt 0 ]; then usage heatmap_reads; exit 1; fi

if [ -f $comp_dir/${names[0]}_${names[1]}_peaks_processed.bed ]; then 
local peaks=$comp_dir/${names[0]}_${names[1]};
else
	if [ -f $comp_dir/${names[0]}_${names[1]}/${names[0]}_${names[1]}_peaks_processed.bed ]; then 	
		local peaks=$comp_dir/${names[0]}_${names[1]}/${names[0]}_${names[1]}
	else
		if [ -f $comp_dir/${names[1]}_${names[0]}_peaks_processed.bed ]; then 
			local peaks=$comp_dir/${names[1]}_${names[0]}
		else
			if [ -f $comp_dir/${names[1]}_${names[0]}/${names[1]}_${names[0]}_peaks_processed.bed ]; then 
				local peaks=$comp_dir/${names[1]}_${names[0]}/${names[1]}_${names[0]}; 
			else
				echo "ERROR: processed file not found, please be sure to specify correct directory for -cp argument"
				usage heatmap_reads;exit 1;
			fi
		fi
	fi
fi
mkdir -p $out_dir

if [ $window_size -gt 0 ]; then 
awk -v len=$window_size -v OFS="\t" '{print $1, int($2+($3-$2)/2-len), int($2+($3-$2)/2+len),$4}' ${peaks}_peaks_processed.bed > $out_dir/${names[0]}_${names[1]}_peaks_processed.bed
else
cp ${peaks}_peaks_processed.bed $out_dir/${names[0]}_${names[1]}_peaks_processed.bed
fi
for ((i=0;i<${#bedgraphs[@]};i++))
do
	$bdg_to_bw ${bedgraphs[i]} $size_file $out_dir/${names[i]}.bw
done
processed=$out_dir/${names[0]}_${names[1]}_peaks_processed.bed
local names=`join_by " " "${names[@]}"`

python $heatmap_mk -b $out_dir -p $processed -r $out_dir -s $order -n $names


}


spacing_impact(){

while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local names=("${!2}")
			echo "names set to: ${2}";shift 2;;
		-cp)
			local comp_dir=${2}
			echo "initial comparison directory set to: ${2}";shift 2;;
		-od)
			local out_dir=${2}
			echo "output directory set to: ${2}";shift 2;;
		-g)
			local genome=${2}
			echo "Genome fasta set to: ${2}";shift 2;;
		-m)
			local tffm=${2}
			echo "matrix set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
		-sth)
			local starting_th=$2
			echo "starting threshold set to: ${2}";shift 2;;
		-eth)
			local ending_th=$2
			echo "ending threshold set to: ${2}";shift 2;;
		-ith)
			local increment_th=$2
			echo "threshold increments set to: ${2}";shift 2;;
		-sp)
			local spacing=("${!2}")
			echo "spacings set to: ${2}";shift 2;;
		-c)
			local colors=("${!2}")
			echo "colors set to: ${2}";shift 2;;
		-h)
			usage spacing_impact; exit;;
		--help)
			usage spacing_impact; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage spacing_impact; exit;;
	esac
done

local Errors=0
if [ -z $names ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $comp_dir ]; then echo "ERROR: -cp argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
if [ -z $tffm ]; then echo "ERROR: -m argument needed";Errors+=1;fi
if [ -z $spacing ]; then echo "ERROR: -sp argument needed";Errors+=1;fi

if [ -z $maxSpace ]; then echo "-maxs argument not used, computing up to 50bp spacing";local maxSpace=50 ;fi
if [ -z $minSpace ]; then echo "-mins argument not used, starting with spacing 0";local minSpace=0 ;fi
if [ -z $offset_left ]; then echo "-maxs argument not used, no left offset used";local offset_left=0 ;fi
if [ -z $offset_right ]; then echo "-maxs argument not used, no right offset used";local offset_right=0 ;fi
if [ -z $colors ]; then echo "-c argument not used, using default palette"; local colors=('#40A5C7' '#307C95' '#205364');fi
# TODO check for error in th
if [ $Errors -gt 0 ]; then usage spacing_impact; exit 1; fi

mkdir -p $out_dir
sed "1d" ${comp_dir}/${names[0]}_${names[1]}/table_${names[0]}_${names[1]}.csv > $out_dir/table_${names[0]}_${names[1]}.bed
local table=$out_dir/table_${names[0]}_${names[1]}.bed

awk -v OFS="\t" '{print $1,$2,$3}' $table > $out_dir/${names[0]}_${names[1]}.bed
local peaks=$out_dir/${names[0]}_${names[1]}.bed

bedtools getfasta -fi $genome -fo $out_dir/${names[0]}_${names[1]}.fa -bed $peaks
local peak_file_fas=$out_dir/${names[0]}_${names[1]}.fa

# python $spacing_mk -tffm $tffm  -o $out_dir/ -n ${names[0]}_${names[1]}_spacing -minInter $minSpace -maxInter $maxSpace  -pos $peak_file_fas  -th 0.01 -neg $peak_file_fas -points True -ol $offset_left -or $offset_right -one_panel  -write_inter -no_absolute_panel
local spacing=`join_by "," "${spacing[@]}"`
echo $spacing
local files=()
local th1=$starting_th
local continu1=1
while [ $continu1 -eq 1 ];
do
	echo -e "th\t$th1"
	sed 's/[:,/]/\t/g' $out_dir/${names[0]}_${names[1]}_spacing.bed | awk -v spacing=$spacing -v th1=$th1 -v th2=$th1 -v OFS="\t" '{save=0; split(spacing,list_spacing,","); for(i=4;i<=NF;i+=3){j=i+1;k=i+2;if(($j>=th1 && $k>=th2)||($j>=th2 && $k>=th1)){gsub(/^ER|^IR|^DR/,"",$i);for(h in list_spacing){if($i == list_spacing[h]){save=1}}}};print $1,$2,$3,save}' | awk -v FS="\t" '{print $4}' > $out_dir/Spacing${th1}.inter

	paste $table $out_dir/Spacing${th1}.inter | sed "1ichr\tbegin\tend\t${names[0]}\t${names[1]}\tSpacing" >  $out_dir/table_${th1}.csv
	files+=("$out_dir/table_${th1}.csv")
	local th1=$(calc $th1+$increment_th)
	local continu1=$(awk -vn1="$th1" -vn2=$ending_th 'BEGIN{print (n1<=n2)?1:0}')
	
done

local files=`join_by "," "${files[@]}"`
local colors=`join_by "," "${colors[@]}"`
Rscript $heatmap_mk -f $files -n "${names[0]},${names[1]}" -d $out_dir -c $colors

}


add_coverage(){

while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-b)
			local bedgraphs=("${!2}")
			echo "bedgraphs to add coverage to the table: ${2}";shift 2;;
        -n)
			local names=("${!2}")
			echo "names of added exp to: ${2}";shift 2;;
        -t)
			local table=$2
			echo "table ouput by initial_comp set to: ${2}";shift 2;;
        -od)
			local out_dir=$2
			echo "name of out directory set to: ${2}";shift 2;;
		-h)
			usage add_coverage; exit;;
		--help)
			usage add_coverage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage add_coverage; exit;;
	esac
done

local Errors=0
if [ -z $bedgraphs ]; then echo "ERROR: -b argument needed";Errors+=1;fi
if [ -z $names ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $table ]; then echo "ERROR: -t argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ $Errors -gt 0 ]; then usage add_coverage; exit 1; fi

awk -v OFS="\t" '{print $1, $2, $3, $4}' $table > $out_dir/tmp_peaks.bed
local peak_file=$out_dir/tmp_peaks.bed
local tmp_table=$out_dir/tmp_table.bed
cp $table $tmp_table

local list_files=()
local i=0
for bdg in ${bedgraphs[@]};
do
	echo $bdg
	echo ${names[$i]}
	bedtools intersect -a $peak_file -b $bdg -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/${names[$i]}.inter
	python $compute_coverage -i $out_dir/${names[$i]}.inter -m
	awk  '{print $5/($3-$2)}' $out_dir/${names[$i]}.inter.cov | sed "1i${names[$i]}" > $out_dir/${names[$i]}.inter.cov.normed
	list_files+=("$out_dir/${names[$i]}.inter.cov.normed")
	local i=$i+1
done
local files=`join_by " " "${list_files[@]}"`
#paste $tmp_table $files > $table
paste $tmp_table $files > $out_dir/table.tsv
rm $peak_file $tmp_table

}


add_score(){

while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-m)
			local matrices=("${!2}")
			echo "gff file of the genome set to: ${2}";shift 2;;
        -n)
			local names=("${!2}")
			echo "osize for the promoters set to: ${2}";shift 2;;
        -t)
			local table=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -od)
			local out_dir=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -od)
			local pocc=True
			echo "Fasta of the genome set to: ${2}";shift 2;;
		-h)
			usage add_coverage; exit;;
		--help)
			usage add_coverage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage add_coverage; exit;;
	esac
done

local Errors=0
if [ -z $matrices ]; then echo "ERROR: -m argument needed";Errors+=1;fi
if [ -z $names ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $table ]; then echo "ERROR: -t argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ $Errors -gt 0 ]; then usage add_coverage; exit 1; fi

awk -v OFS="\t" '{print $1, $2, $3, $4}' $table > $out_dir/tmp_peaks.bed
local peak_file=$out_dir/tmp_peaks.bed
local tmp_table=$out_dir/tmp_table.bed
cp $table $tmp_table

local list_files=()
local i=0
for matrix in $matrices[@];
do
	if [[ $matrice == *".xml"* ]]; then # TFFM scores computation
		python $tffmscores -o ${out_dir}/tffm_scores.tsv -pos $peak_fasta -t ${matrice}
		awk '{print $8}' ${out_dir}/${name[$i]}_scores.tsv 
		list_files+=("${out_dir}/${name[$i]}_scores.tsv")
	fi
	if [[ $matrice == *".txt"* ]]; then # K-mer scores computation
		echo "WIP"
	fi
	if [[ $matrice == *".pfm"* ]]; then  # PWM scores
		python $scores_prog -m ${matrice} -f $peak_fasta -o ${out_dir}/
		sort -u -k1,1 -k4,4nr ${out_dir}/$(basename $peak_fasta).scores | sort -u -k1,1 | awk 	'{print $4}' | sort -nr >"${out_dir}/tab_${name[$i]}.tsv"
		list_files+=("${out_dir}/tab_${name[$i]}.tsv")
		if $pocc ; then # PWM Pocc
			python $pocc_pfm -s ${out_dir}/$(basename $peak_fasta).scores -o ${out_dir}/Pocc.pocc
			sort -nr ${out_dir}/Pocc.pocc > ${out_dir}/${name[$i]}_Pocc_scores.tsv
			list_files+=("${out_dir}${name[$i]}_Pocc_scores.tsv")
		fi
	fi
	local i=${i}+1
done
local files=`join_by " " "${list_files[@]}"`
paste $tmp_table $files > $table

rm $peak_file $tmp_table

}

