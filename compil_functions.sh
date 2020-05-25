# list passed as arg
# list=("a" "b")
# test(){ testVar=("${!1}"); echo "${testVar[@]}";  }
# test list[@]

source /home/304.6-RDF/scripts/DAP_global_analysis_p3.7/compil_usages.sh


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


# WIP
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
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage mapping_Fastq_bowtie2; exit;;
	esac
done

if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; local proc=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; usage; exit 1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; usage; exit 1; fi

PATH_TO_BOWTIE_INDEX="/home/304.3-STRUCTPDEV/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"
PATH_TO_BOWTIE2=$(dirname $(which bowtie2))
PATH_TO_SAMTOOLS=$(dirname $(which samtools))
PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools)))

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

# WIP
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
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage main_mapping_Fastq; exit;;
	esac
done

if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; local proc=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $Fastq2 ]; then echo "-f2 argument not used, assuming non paired-end analysis"; local Fastq2=""; fi
if [ -z $Fastq_dir ]; then echo "ERROR: -fd argument needed"; usage; exit 1; fi
if [ -z $Mapping_dir ]; then echo "ERROR: -md argument needed"; usage; exit 1; fi
if [ -z $Fastq1 ]; then echo "ERROR: -f1 argument needed"; usage; exit 1; fi
if [ -z $Name ]; then echo "ERROR: -n argument needed"; usage; exit 1; fi

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


# WIP
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
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage peakcalling_MACS2; exit;;
	esac
done

if [ -z $Genome_length ]; then echo "-g argument not used, assuming A.thaliana is used: 120000000"; local proc=120000000; fi
if [ -z "$additionnal_args" ]; then echo "-add argument not used, no additional argument passed to MACS2"; local additionnal_args=""; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; usage; exit 1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; usage; exit 1; fi


PATH_TO_MACS=$(dirname $(readlink -f $(which macs2)))
PATH_TO_SAMTOOLS=$(dirname $(which samtools))
PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools)))

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
		-h)
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage main_peakcalling; exit;;
	esac
done

if [ -z $Genome_length ]; then echo "-g argument not used, assuming A.thaliana is used: 120000000"; local proc=120000000; fi
if [ -z "$additionnal_args" ]; then echo "-add argument not used, no additional argument passed to MACS2"; local additionnal_args=""; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; usage; exit 1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; usage; exit 1; fi
if [ -z $input_dir ]; then echo "ERROR: -cd argument needed"; usage; exit 1; fi
if [ -z $name_cons ]; then echo "ERROR: -nc argument needed"; usage; exit 1; fi

local list_peaks_mspc=()
local list_bdg=()

for dir in ${in_dir[@]};
do
  for bam in $(find $input_dir -name "*.filtered.sorted.nodup.bam")
  do
    cp $bam $dir/control.filtered.sorted.nodup.bam
  done
  local name=${dir##*/}
  list_peaks_mspc+=("-i $out_dir/$name/${name}_peaks.bed")
  list_bdg+=("$out_dir/$name/${name}_cpm.bdg")
#   peakcalling_MACS2 $dir $out_dir $Genome_length "$additionnal_args"
   peakcalling_MACS2 -id $dir -od $out_dir -g $Genome_length -add "$additionnal_args"
done
local list_peaks_mspc=`join_by " " "${list_peaks_mspc[@]}"`

if [[ $(uname -a) =~ "el7" ]]; then
  mspc ${list_peaks_mspc} -r Tec -w 1e-4 -s 1e-8 -c 100% -o $out_dir/$name_cons
else
  echo "Error SL7 node is needed for mspc" && exit 1
fi

sed '1d' $out_dir/$name_cons/ConsensusPeaks.bed | awk -v OFS="\t" '{print $1,$2,$3}' > $out_dir/$name_cons/$name_cons.bed

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
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage; exit;;
	esac
done

if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; usage; exit 1; fi
if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; usage; exit 1; fi
if [ -z $result ]; then echo "ERROR: -od argument needed"; usage; exit 1; fi
if [ -z $data ]; then echo "ERROR: -id argument needed"; usage; exit 1; fi
if [ -z $filterCov ]; then echo "no filter on coverage applied"; local filterCov=0; fi
if [ -z $filterHeight ]; then echo "no filter on peak height applied"; local filterHeight=0; fi

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

	export TMPDIR=/nobackup
	local do_plot_quant=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/Hist_cov_gen.r
	local merge_peaks=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/merge_all_peaks.py
	local compute_coverage=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/compute_coverage.py

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
			echo "data directory set to: ${2}";shift 2;;
		-n)
			local name=$2
			echo "data directory set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "data directory set to: ${2}";shift 2;;
		-od)
			local results_motif=$2
			echo "data directory set to: ${2}";shift 2;;
		-ls)
			local learning_size=$2
			echo "data directory set to: ${2}";shift 2;;
		-nm)
			local nb_motifs=$2
			echo "data directory set to: ${2}";shift 2;;
		-mt)
			local motifstffm=$2
			echo "data directory set to: ${2}";shift 2;;
		-mim)
			local min_motif=$2
			echo "data directory set to: ${2}";shift 2;;
		-mam)
			local max_motif=$2
			echo "data directory set to: ${2}";shift 2;;
		-pal)
			local pal=true
			echo "palindromic mode activated";shift 1;;
		-s)
			local seed=$2
			echo "data directory set to: ${2}";shift 2;;
		-h)
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage; exit;;
	esac
done

if [ -z $peaks ]; then echo "ERROR: -p argument needed"; usage; exit 1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; usage; exit 1; fi
if [ -z $results_motif ]; then echo "ERROR: -od argument needed"; usage; exit 1; fi
if [ -z $learning_size ]; then echo "-ls argument not used, using 600 peaks"; local  learning_size=600; fi
if [ -z $motifstffm ]; then echo "-mt argument not used, first motif used for tffm creation"; local motifstffm=0; fi
if [ -z $min_motif ]; then echo "-mim argument not used, min size of motif is 6"; local min_motif=6; fi
if [ -z $max_motif ]; then echo "-mam argument not used, max size of motif is 15"; local min_motif=15; fi
if [ -z $nb_motifs ]; then echo "-nm argument not used, searching for 1 motif"; local nb_motifs=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if ! $pal ; then echo "palindromic mode not activated"; fi
if [ $max_motif -lt $min_motif ]; then
	echo "WARNING: things got mixed up, max size of motif can't be lower than min size."
	echo -e "\tUsing min size as max size "
	local max_motif=$(calc $max_motif+$min_motif) # max motif = total
	local min_motif=$(calc $max_motif-$min_motif) # new min_motif = old max_motif (total - old min = old max)
	local max_motif=$(calc $max_motif-$min_motif) # new max_motif = old min_motif (total - new min = old min)
fi

meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip
meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme
meme2pfm=/home/304.6-RDF/scripts/meme2pfm.sh
# meme2pfm=/home/304.6-RDF/scripts/DAP_global_analysis/meme2pfm_modifJM.sh
prepMEMEforPalTFFM=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/prepMEMEforPalTFFM.py
pfmTOtffm=/home/304.6-RDF/LFY/scripts/get_tffm.py

mkdir -p $results_motif/$name/sets $results_motif/$name/meme $results_motif/$name/tffm
    
    head -$learning_size $peaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
    bedtools getfasta -fi $genome -fo $results_motif/$name/sets/${name}.fas -bed $results_motif/$name/sets/${name}.bed
    
    if $pal ; then
    echo "pal"
    $meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-pal -meme-maxsize $(calc $learning_size*1000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs 1 -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed 
    else
    $meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-maxsize $(calc $learning_size*1000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs 1 -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed
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
			echo "output directory set to: ${2}";shift 2;;
        -p)
			local promoterLength=$2
			echo "output directory set to: ${2}";shift 2;;
        -g)
			local Genome=$2
			echo "output directory set to: ${2}";shift 2;;
		-h)
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage; exit;;
	esac
done

prepare_gff=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/prepare_gff.py
generate_BedFromGff=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/generate_bed_gff.py

python $prepare_gff -gi ${bed_genome}.gff3 -go ${bed_genome}_prep.gff3  -f ${Genome}
python $generate_BedFromGff -g ${bed_genome}_prep.gff3 -p $promoterLength -o ${bed_genome}_prep.bed

bedtools sort -i ${bed_genome}_prep.bed | bedops --partition - | bedmap --echo --echo-map-id --delim '\t' - ${bed_genome}_prep.bed | awk -v OFS="\t" '$4==""{print $1,$2,$3,"intergenic";next} {print $0}' > ${bed_genome}.bed


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
			echo "data directory set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "data directory set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-anf)
			local annotation_file=$2
			echo "output directory set to: ${2}";shift 2;;
		-nb)
			local number_of_NS=$2
			echo "output directory set to: ${2}";shift 2;;
		-bs)
			local bin_size=$2
			echo "output directory set to: ${2}";shift 2;;
		-gc)
			local deltaGC=$2
			echo "output directory set to: ${2}";shift 2;;
		-lt)
			local limit_type=$2
			echo "output directory set to: ${2}";shift 2;;
		-s)
			local seed=$2
			echo "output directory set to: ${2}";shift 2;;
		-h)
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage; exit;;
	esac
done

if [ -z $peaks ]; then echo "ERROR: -p argument needed"; usage; exit 1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; usage; exit 1; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/304.6-RDF/data/tair10.fas"; local genome="/home/304.6-RDF/data/tair10.fas" ; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; usage; exit 1; fi
if [ -z $annotation_file ]; then echo "ERROR: -anf argument needed"; usage; exit 1; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, searching for 1 negative set"; local number_of_NS=1; fi
if [ -z $bin_size ]; then echo "-bs argument not used, using 250bp as bin size"; local bin_size=250 ; fi
if [ -z $deltaGC ]; then echo "-gc argument not used, using 0.03 as maximum of GC% divergence"; local deltaGC=0.03 ; fi
if [ -z $limit_type ]; then echo "-lt argument not used, using 1000 as minimum group size for origin (use -h for more information)"; local limit_type=1000 ; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi

negative_set_script=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/negative_set_generator_multiThread.py

mkdir -p ${results}
# awk -v OFS='\t' ' ($3-$2)>1200 {print $1, ($2+($3-$2)-600), ($2+($3-$2)+600);next} {print $1, $2, $3}' $peak > ${results}/${name}_fil.bed
python $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $bin_size

}


compute_ROCS(){

# compute_ROCS -p -ns -m -n -g -od -pc
local pocc=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks file set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "peaks file set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "peaks file set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "data directory set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "data directory set to: ${2}";shift 2;;
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
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage; exit;;
	esac
done

pocc_pfm=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/compute_POcc.py
tffmscores=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/get_best_score_tffm.py
scores_prog=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/scores.py
plot_ROCS_prog=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/plots_ROCS_multiple.py

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
			echo "peaks file set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "peaks file set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "peaks file set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "data directory set to: ${2}";shift 2;;
        -th)
			local thresholds=("${!2}")
			echo "data directory set to: ${2}";shift 2;;
        -maxy)
			local maxy=$2
			echo "data directory set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "data directory set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "data directory set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "data directory set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "data directory set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "data directory set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "data directory set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-pc)
			local pocc=true
			echo "Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
		-h)
			usage; exit;;
		--help)
			usage; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage; exit;;
	esac
done
mkdir -p $results
spacing_mk=/home/304.6-RDF/scripts/DAP_global_analysis_p3.7/get_interdistances.py

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




