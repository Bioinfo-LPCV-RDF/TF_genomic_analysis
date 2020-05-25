

usage(){
	case $1 in
	mapping_Fastq_bowtie2)
echo -e "
==========
usage: mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s [int] -pr [int]

general infos:

-- Mandatory arguments:
    -id     PATH    :    Set the directory where the Fastq is stored
                         In case of paired-end analysis, please regroup your 
                         Fastq in the same directory.
    -od     PATH    :    Set the Output directory, where you will find the BAM
                         files created by this programs.
    
-- optional arguments :
    -s      INT     :    Set the seed for the random (default: 1254)
    -pr     INT     :    Set the number of threads to use. (default: 1)
\n";;
	main_mapping_Fastq)
echo -e "
==========
usage: main_mapping_Fastq -fd <PATH> -md <PATH> -f1 <FILE> -f2 [FILE] 
       -n <STRING> -s [INT] -pr [INT]

general infos:

-- Mandatory arguments:
    -fd     PATH    :    Set the directory where the Fastq is stored. If the
                         Fastq are stored in another folder, they will be copy
                         pasted in the directory selected here.
    -md     PATH    :    Set the Output directory, where you will find the BAM
                         files created by this programs.
    -f1     FILE    :    Set the Fastq file to be used for the mapping. In case
                         of paired-end analysis please separate pairs in two
                         files (see -f2 argument)
    -n      STRING  :    Name of the sub-folder to be created in the output
                         directory.

-- optional arguments :
    -f2     FILE    :    Set the Fastq file for the second mate of the pairs.
                         leave empty in case of single end analysis.
    -s      INT     :    Set the seed for the random (default: 1254)
    -pr     INT     :    Set the number of threads to use. (default: 1)
\n";;
	peakcalling_MACS2)
echo -e "
==========
usage: peakcalling_MACS2 -id <PATH> -od <PATH> -g <FILE> -add [STRING]

general infos:

-- Mandatory arguments:
    -id     PATH    :    Set the directory where the BAM files are stored
    -od     PATH    :    Set the Output directory, where you will find the 
                         BEDGRAPH and NarrowPeak files created by this 
                         programs.

-- optional arguments :
    -g      INT     :    Set the seed for the random (default: 1254)
    -add     INT     :    Set the number of threads to use. (default: 1)
\n";;
	main_peakcalling) # main_peakcalling -id -cd -od -nc -g -add
echo -e "
==========
usage: main_peakcalling -id <PATH> -cd <PATH> -od <PATH> -nc <STRING> -g [FILE]
       -add [STRING]

general infos:

-- Mandatory arguments:
    -id     PATH    :    Set the directory where the BAM files are stored
    -cd     PATH    :    Set the directory where the BAM files for the control
                         sets are stored
    -od     PATH    :    Set the Output directory, where you will find the 
                         BEDGRAPH and NarrowPeak files created by this 
                         programs.
    -nc     STRING  :    Set the name of the final output directory where the
                         consensus NarrowPeak & BEDGRAPH will be stored

-- optional arguments :
    -g      INT     :    Set the seed for the random (default: 1254)
    -add     INT     :    Set the number of threads to use. (default: 1)
\n";;
	initial_comparison) # main_peakcalling -id -cd -od -nc -g -add
echo -e "
==========
usage: initial_comparison -n1 <STRING> -n2 <STRING> -od <PATH> -id <PATH> 
       -r1 <INT> -r2 <INT> -f [INT]

general infos:

-- Mandatory arguments:
    -n1     STRING  :    name of directory for dataset 1
    -n2     STRING  :    name of directory for dataset 2
    -od     PATH    :    Set the Output directory, where you will find the 
                         BEDGRAPH and NarrowPeak files created by this 
                         programs.
    -id     PATH    :    Set the directory where the BED/BEDGRAPH files are
                         stored.
    -r1     INT     :    Set the ratio of CFR from which the peaks will be 
                         considered specific of the datatset 1
    -r2     INT     :    Set the ratio of CFR from which the peaks will be 
                         considered specific of the datatset 2

-- optional arguments :
    -f      INT     :    Set the threshold of coverage for both sample in 
                         order for the peaks to be considered as significant
\n";;
	compute_motif) 
echo -e "
==========
usage: compute_motif -p <FILE> -n <STRING> -od <PATH> -g [FILE] -ls [INT]
       -nm [INT] -mim [INT] -mam [INT] -mt [INT] [-pal] -s [INT]

general infos:

-- Mandatory arguments:
    -p     FILE     :    name of directory for dataset 1
    -n     STRING   :    name of directory for dataset 2
    -od    PATH     :    Set the Output directory

-- optional arguments :
    -g      FILE    :    Set the genome file to use as reference
    -ls     INT     :    Set the learning-size for both PWM and TFFM 
                         matrices
    -nm     INT     :    Set the number of PWM matrices to compute
    -mim    INT     :    Set the minimum width of PWM matrices
    -mam    INT     :    Set the maximum width of PWM matrices
    -mt     INT     :    Set which PWM matrix to use as model for TFFM 
                         computation (starts as 0)
    -pal            :    Activate the palindromic mode for PWM matrices. 
                         Does not work for TFFM computation
    -s      INT     :    Set the seed for random
\n";;
	prep_annotation) 
echo -e "
==========
usage: prep_annotation -n <FILE> -p <INT> -g <FILE>

general infos:

-- Mandatory arguments:
    -n     FILE     :    gff file of the genome, stripped of its extension
    -p     STRING   :    size for the promoters
    -g     FILE     :    Fasta of the genome used as reference
\n";;
	compute_NS) 
echo -e "
==========
usage: compute_NS -p <FILE> -n <STRING> -g <FILE> -od <PATH> -anf <FILE> -nb [INT] -bs [INT] -gc [INT] -lt [INT] -s [INT]

general infos:

-- Mandatory arguments:
    -n     STRING   :    name for the results files
    -p     FILE     :    peak file (bed format)
    -g     FILE     :    Fasta of the genome used as reference
    -od    PATH     :    Set the Output directory
    -anf   PATH     :    Annotation file, created by prep_annotation function

-- optional arguments :
    -nb    INT      :    number of separate negative sets file to created
    -bs    INT      :    size of the bin used. Ideally the  bin size should
                         be set 30-50 higher than the smallest peak in the 
                         positive set
    -gc    INT      :    ideal maximum difference wanted for %GC. Note that 
                         it can be moved higher by the program if need be (not
                         enough regions with similar GC%)
    -lt    INT      :    minimum number of genomic regions sharing the same 
                         complex origin (example : intron;exon). 
    -s     INT      :    Set the seed for random
\n";;
	compute_ROCS) 
echo -e "
==========
usage: compute_ROCS -p <FILE1> <FILEn> -ns <FILE1> <FILEn> -m <FILE1> <FILEn>
       -n <STRING1> <STRINGn> -od <PATH> -g [FILE] [-pc]

general infos:

-- Mandatory arguments:
    -ns    FILEs    :    name of the negative sets
    -p     FILEs    :    peak file (bed format)
    -m     FILEs    :    matrices files
    -od    PATH     :    Set the Output directory
    -n     STRING   :    name for the sets in the ROCs

-- optional arguments :
    -g     FILE     :    Fasta of the genome used as reference
    -pc             :    pocc mode. 
\n";;
	compute_space) 
echo -e "
==========
usage: compute_space -p <FILE1> <FILEn> -ns <FILE1> <FILEn> -m <FILE1> <FILEn>
       -n <STRING1> <STRINGn> -od <RESULT_DIR> -th <INT1> <INTn> -maxy [INT]
       -maxs [INT] -mins [INT] -ol [INT] -or [INT] -g [FILE] -nb <NB_of_NS>
       [-pc]

general infos:

-- Mandatory arguments:
    -p     FILEs    :    peak file (bed format)
    -ns    FILEs    :    name of the negative sets
    -m     FILEs    :    matrices files
    -n     STRINGs  :    name for the sets in the ROCs
    -od    PATH     :    Set the Output directory
    -th    FLOATs   :    list of thresholds to use. (between 0-1 for TFFM
                         matrices & between -60-0 for PWMs

-- optional arguments :
    -maxy  INT      :    
    -maxs  INT      :    
    -mins  INT      :    
    -ol    INT      :    
    -or    INT      :    
    -g     FILE     :    
    -pc             :    
\n";;
	*)
echo -e "
usage: call one of the function below to see more details about this specific function

--functions implemented:
mapping_Fastq_bowtie2 :  
main_mapping_Fastq    :  
peakcalling_MACS2     :  
main_peakcalling      :  
initial_comparison    :  
compute_motif         :  
prep_annotation       :  
compute_NS            :  
compute_ROCS          :  
compute_space         :  
\n";;
	esac
echo -e "==========\n"
}
















