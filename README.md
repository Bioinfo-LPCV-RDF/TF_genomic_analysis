  
This pipeline is used for DAP-seq analysis of MADS transcription Factors.  
Compil_function.sh gathers all the functions needed for the analysis and Wrapper_MADS_tetra_paper.sh call all these functions, with specific parameters.   
The codes to generate the figures 5, S8 and S10 are not in the wrapper and have to be launched separately (see ### 10. and ### 11. ).  
  
## DEPENDENCIES  
  
<pre>  
fastqc #v0.11.7 mapping_Fastq_bowtie2  
NGmerge #v0.2_dev mapping_Fastq_bowtie2  
bedtools #v2.27.1 peakcalling_MACS2  
mspc #v4.0.0 main_peakcalling  
inkscape #v0.92.2 compute_space  
python v2.7.9  #Python 3 version in progress  
    argparse  
	tffm_module  
	Bio.Seq  
	pybedtools  
	numpy  
	tempfile  
	matplotlib  
	sklearn.metrics  
	metaseq  
	multiprocessing  
R v3.5.0  
	Biostrings  
	ggplot2  
	Cairo  
	latex2exp  
	optparse  
</pre>  
  
## INSTALLATION  
<pre>  
please be sure to install all dependencies  
  
be sure to have repositories structured as below:  
  
/  
compil_functions.sh  
compil_usages.sh  
bin/  
	compute_coverage.py  
	compute_Pocc.py  
	compute_spacingScorev2.r  
	generate_bed_gff.py  
	get_best_score_tffm.py  
	get_interdistances.py  
	get_tffm.py  
	Hist_cov_gen.r  
	interdistances_functions.py  
	meme2pfm.sh  
	merge_all_peaks.py  
	negative_set_generator.py  
	plot_functions.py  
	plots_ROCS_multiple.py  
	prepare_gff.py  
	prepMEMEforPalTFFM.py  
	scores.py  
	show_reads.py  
  
if you work on arabidopsis, the annotation file needed for negative sets generation is provided (named tair10.bed)  
You can change the pathways in the beginning of <compil_functions.sh>  
</pre>  
## PIPELINE  
  
### 1.Mapping  
  
Using bowtie2 and single-end or paired-end Fastq.gz.  
  
usage:   
> main_mapping_Fastq -fd [PATH] -md [PATH} -f1 [FILE} -f2 [FILE]   
>       -n [STRING} -s [INT] -pr [INT]  
>  
general infos: directory and file preparation for mapping step  
>  
-- Mandatory arguments:  
>    -fd     PATH    :    Set the directory where the Fastq is stored. If the  
>                         Fastq are stored in another folder, they will be copied  
>                         in the directory selected here.  
>    -md     PATH    :    Set the Output directory, where you will find the BAM  
>                         files created by this programs.  
>    -f1     FILE    :    Set the Fastq file to be used for the mapping. In case  
>                         of paired-end analysis, please separate the pairs in two  
>                         files (see -f2 argument)  
>    -n      STRING  :    Name of the sub-folder to be created in the output  
>                         directory.  
>  
-- optional arguments :  
>    -f2     FILE    :    Set the Fastq file for the second mate of the pairs.  
>                         leave empty in case of single end analysis.  
>    -s      INT     :    Set the random seed (default: 1254)  
>    -pr     INT     :    Set the number of threads to use. (default: 1)  
>     
### 2.Peak Calling  
  
Using Macs2 and MSPC to treat multiple replicates. After MSPC, peaks are resized to 400pb around maximums if they are separated by more than 200pb. This is to ensure that two close peaks are not considered as one.  
  
example for peak resizing:  
</>: extremities for the peaks  
| : maximum of peaks, determined by Macs2  
position                          200                       605         
replicate 1:         <-------------|------------------->  
replicate 2              <---------|----------------------->  
replica:                            <------------------------|-------->  
  
consensus MSPC:      <------------------------------------------------>  
peaks resized          <-----------|-----------> <-----------|----------->  
  
  
usage:   
> main_peakcalling -id [PATH} -cd [PATH} -od [PATH} -nc [STRING}   
>       -g [FILE] -add [STRING] -top [INT]  
>  
>general infos: file/directory preparation, peakcalling_MACS2 launching,  
>               Multiple Sample Peak Calling (MSPC) procedure, union of   
>               coverages, peak resizing.  
>  
-- Mandatory arguments:  
>    -id     PATH    :    Set the directory where the BAM files are stored  
>    -cd     PATH    :    Set the directory where the BAM files for the control  
>                         sets are stored  
>    -od     PATH    :    Set the Output directory, where you will find the   
>                         BEDGRAPH and NarrowPeak files created by this   
>                         programs.  
>    -nc     STRING  :    Set the name of the final output directory where the  
>                         consensus NarrowPeak & BEDGRAPH will be stored  
>  
-- optional arguments :  
>    -g      INT     :    Set the genome length to use as reference for peak   
>                         calling  
>    -add    STRING  :    Set the additionnal parameters for MACS2 (option used  
>                         by default: -t -c -B -f -n -g; please do not   
>                         overwrite those)  
>    -top    INT     :    Set the maximum number of peaks to use. Peaks will be  
>                         ordered by p-value before the cut-off is applied  
  
### 3.COMPARISONS of samples/replicates  
  
  
  
usage:   
> initial_comparison -n1 [STRING] -n2 [STRING] -od [PATH] -id [PATH]   
>       -r1 [INT] -r2 [INT] -f [INT]  
>  
<pre>  
general infos: comparison of peak coverages between two samples or replicates,   
               filtered (optional) by the max or the coverage of each peak,   
               plot the coverage of each peak in the first sample against its   
               coverage in the other sample.   
               peak reattribution according to CFR (and not according to peak caller;   
               plots with both for comparison, CFR: coverages fold reduction)  
</pre>  
  
-- Mandatory arguments:  
>    -n1     STRING  :    name of directory for dataset 1  
>    -n2     STRING  :    name of directory for dataset 2  
>    -od     PATH    :    Set the Output directory, where you will find the   
>                         BEDGRAPH and NarrowPeak files created by this   
>                         programs.  
>    -id     PATH    :    Set the directory where the BED/BEDGRAPH files are  
>                         stored.  
>    -r1     INT     :    Set the coverage fold ratio -CFR- from which the   
>                         peaks will be considered specific of the datatset 1  
>    -r2     INT     :    Set the CFR from which the peaks will be   
>                         considered specific of the datatset 2  
>  
-- optional arguments :  
>    -f      INT     :    Threshold on coverage -in RPKM-. Filters out peaks  
>                         for which for both samples do not pass the threshold.  
>    -he     INT     :    Threshold on the maximum height of coverage -in RPKM-  
>                         Filters out peaks that have a max   
>                          coverage below this value in both samples.  
  
### 4.MOTIF DISCOVERY  
  
  
  
usage:   
>compute_motif -p [FILE} -n [STRING] -od [PATH] -g [FILE] -ls [INT]  
>       -nm [INT] -mim [INT] -mam [INT] -mt [INT] [-pal] -s [INT]  
>  
<pre>  
general infos: PWM and TFFM motif generation using learning set (creation of   
              testing set for later use).  
			  - PWM: Position Weight Matrix  
			  - TFFM: Transcription Factor Flexible Model  
</pre>  
>  
-- Mandatory arguments:  
>    -p     FILE     :    peak file to use (bed, narrowpeak format accepted)  
>    -n     STRING   :    name for result directory (to be created)  
>    -od    PATH     :    Set the Output directory (will include the   
>                         subdirectory created by -n option)  
>  
-- optional arguments :  
>    -g      FILE    :    Set the genome file to use as reference  
>    -ls     INT     :    Set the learning-size for both PWM and TFFM   
>                         matrices  
>    -nm     INT     :    Set the number of PWM matrices to compute  
>    -mim    INT     :    Set the minimum width of PWM matrices  
>    -mam    INT     :    Set the maximum width of PWM matrices  
>    -mt     INT     :    Set which PWM matrix to use as model for TFFM   
>                         computation (starts as 0)  
>    -pal            :    Activate the palindromic mode for PWM matrices.   
>                         Does not work for TFFM computation  
>    -s      INT     :    Set the seed for random  
  
### 5.NEGATIVE SET GENERATION  
  
  
  
>usage: compute_NS -p [FILE] -n [STRING] -g [FILE] -od [PATH] -anf [FILE]   
>       -nb [INT] -bs [INT] -gc [INT] -lt [INT] -s [INT]  
>  
<pre>  
general infos: creation of negative sets (for ROCs and spacing analysis).   
               Negative sequences are created as close as possible to   
               the \"positive\" sequences (the peaks), For each positive sequence,  
               a genomic sequence having the same origin (promoters, exon,   
               intergenic, etc), a close GC percentage, and being not bound  
               in the positive set is generated.   
</pre>  
-- Mandatory arguments:  
>    -n     STRING   :    name for the results files  
>    -p     FILE     :    peak file (bed format)  
>    -g     FILE     :    reference genome fasta file  
>    -od    PATH     :    Set the Output directory  
>    -anf   PATH     :    Annotation file, created by prep_annotation function  
>  
-- optional arguments :  
>    -nb    INT      :    number of separate negative sets file to created   
>                         (default: 1)  
>    -ws    INT      :    size of the window used to slice the genome. Ideally   
>                         the it should be set 30-50 higher than the smallest  
>                         peak in the positive set (default: 250)  
>    -gc    INT      :    ideal maximum difference wanted for %GC. Note that   
>                         it can be moved higher by the program if need be (not  
>                         enough regions with similar GC%; default: 0.03)  
>    -lt    INT      :    minimum number of genomic regions sharing the same   
>                         complex origin (example : intron;exon; default: 1000)  
>    -s     INT      :    Set the seed for random  
  
### 6.ROCS  
  
usage:   
> $ compute_ROCS -p [FILE1] [FILEn] -ns [FILE1] [FILEn] -m [FILE1] [FILEn]  
 >      -n [STRING1] [STRINGn] -od [PATH] -g [FILE] [-pc] -color [STRING1]  [STRINGn]  
  
<pre>  
general infos: computes scores for PWM/TFFM on negative sequences and peaks.  
               Computes ROCs and AUCROC from those scores.  
</pre>  
-- Mandatory arguments:  
>    -ns    FILEs    :    name of the negative sets  
>    -p     FILEs    :    peak file (bed format)  
>    -m     FILEs    :    matrices files  
>   -od    PATH     :    Set the Output directory </pre>  
>    -n     STRING   :    name for the sets in the ROCs  
>  
-- optional arguments :  
>    -g     FILE     :    Fasta of the genome used as reference.  
>    -pc             :    pocc mode.   
>    -color STRINGs  :    list of colors to use in each set (hexadecimal).  
  
### 7.SPACING  
  
usage:  
> $ compute_space -p [FILE1] [FILEn] -ns [FILE1] [FILEn] -m [FILE1] [FILEn]  
>       -n [STRING1] [STRINGn] -od [RESULT_DIR] -th [INT1] [INTn]   
>       -nb {NB_of_NS] -maxy [INT] -maxs [INT] -mins [INT] -ol [INT] -or [INT]  
>       -g [FILE] [-pc]  
  
<pre>  
general infos: computes spacing analysis on peak files and negatives sets   
               using PWM/TFFM.  
</pre>  
-- Mandatory arguments:  
>    -p     FILEs    :    peak files (bed, narrowpeak format)  
>    -ns    FILEs    :    name of the negative sets  
>    -nb    INT      :    number of negatives sets to use against each peak   
>                         file (suggested: 3; be sure to create enough negative  
>                         file using compute_NS)  
>    -m     FILEs    :    matrices files  
>    -n     STRINGs  :    name for the sets in the ROCs  
>    -od    PATH     :    Set the Output directory  
>    -th    FLOATs   :    list of thresholds to use. (between 0-1 for TFFM  
>                         matrices & between -60-0 for PWMs). Choosing   
>                         thresholds can be quite tough. One way to get  
>                         satisfying threshold is to set them empirically and  
>                         adjust them according to the results. If dots are   
>                         missing for some spacing, thresholds are too high,   
>                         reduce them (closer to 0 for TFFM; closer to -60 for   
>                         PWM). If the enrichment increases and reduces itself   
>                         rapidly, everywhere in the graph, you are seeing   
>                         noises, you must increase you thresholds (closer to 1 for   
>                         TFFM; closer to 0 for PWM).  
>  
-- optional arguments :  
>    -maxy  INT      :    maximum enrichment to display on y-axis.(default: NA)  
>    -maxs  INT      :    maximum spacing to compute. For dimers (using   
>                         monomeric matrices), a value between 30 to 50 bp   
>                         seems a good start. For tetramers (using dimeric   
>                         matrices) a value between 70 to 100 bp is a good   
>                         start. You may adjust those values after your first  
>                         results. (default: )  
>    -mins  INT      :    minimum spacing to compute. Usually set to 0. You may  
>                         change this value if you add offsets, this will  
>                         prevent huge drop off in enrichment for the first   
>                         spacings computed.  
>    -ol    INT      :    offset to apply on the left of the matrix used. This  
>                         argument (in addition to -or) allows to change the   
>                         way spacings are counted. This may be usefull when   
>                         you want to count space around a consensus sequence  
>                         or from the center of the matrix.  
>    -or    INT      :    offset to apply on the right of the matrix used. See  
>                         -ol option for details.  
>    -g     FILE     :    Fasta of the genome used as reference.  
>    -pc             :    pocc mode.  
  
### 8.Heatmap  
  
  
  
usage:  
> $ heatmap_reads -b [FILE1] [FILE2] -n [STRING1] [STRING2] -od [PATH]   
>       -cp [PATH] -s [FILE] -or [INT] -ws [INT]   
>  
<pre>  
general infos: computes heatmap of coverages to compare two sets of given peaks.   
               \"initial_comp\" function must have been used before and its   
               results directory specified for \"-cp\" argument.  
</pre>  
-- Mandatory arguments:  
>    -b     FILEs    :      
>    -n     STRINGs  :    name for the sets in the plot  
>    -cp    PATH     :    directory \"initial_comp\" of the two samples  
>    -s     FILE     :    chromosome sizes file (name\\tsize format)  
>    -od    PATH     :    Set the Output directory  
>  
-- optional arguments :  
>    -or     INT     :    choose which variable to use for ordering the   
>                         sequences. 3 choices: 1 to order by 1st set, 2 to   
>                         order by 2nd set and 3 to order by ratio 1st/2nd.  
>                         (default: 3)  
>    -ws     INT     :    choose the size of the window to plot. (default:   
>                         sequence length)  
  
### 9.IMPACT OF SPACING  
  
  
  
usage:   
> $ spacing_impact -n [STRING1] [STRING2] -cp [PATH] -od [PATH] -g [FILE]   
>       -m [FILE] -maxs [INT] -mins [INT] -ol [INT] -or [INT] -sth [INT] -eth   
>       [INT] -ith [INT] -sp [FLOAT1] .. [FLOATN] -c [STRING1] .. [STRINGN]  
>  
<pre>  
general infos: creates a plot showing the relation of the ratio of coverages   
               and the presence of some specific spacing.  
</pre>  
-- Mandatory arguments:  
>    -n     STRINGs  :    name for the sets in the plot  
>    -cp    PATH     :    directory \"initial_comp\" of the two samples  
>    -od    PATH     :    Set the Output directory  
>    -g     FILE     :    Fasta of the genome used as reference.  
>    -m     FILE     :    matrice file  
>    -sp    FLOATs   :    Spacing to search specifically.  
>  
-- optional arguments :  
>    -maxs  INT      :    maximum spacing to compute. For dimers (using   
>                         monomeric matrices), a value between 30 to 50 bp   
>                         seems a good start. For tetramers (using dimeric   
>                         matrices) a value between 70 to 100 bp is a good   
>                         start. You may adjust those values after your first  
>                         results. (default: )  
>    -mins  INT      :    minimum spacing to compute. Usually set to 0. You may  
>                         change this value if you add offsets, this will  
>                         prevent huge drop off in enrichment for the first   
>                         spacings computed.  
>    -ol    INT      :    offset to apply on the left of the matrix used. This  
>                         argument (in addition to -or) allows to change the   
>                         way spacings are counted. This may be usefull when   
>                         you want to count space around a consensus sequence  
>                         or from the center of the matrix.  
>    -or    INT      :    offset to apply on the right of the matrix used. See  
>                         -ol option for details.  
>    -sth   INT      :    lower value for threshold  
>    -eth   INT      :    upper value for threshold  
>    -ith   INT      :    incremental value for threshold  
>    -c     STRINGs  :    sets of color to use in plot (hexadecimal format)  
  
### 10 Crossing DAP-seq, ChIP-seq and gene expression (figure 5 and figure S10)  
  
<pre>   
general infos : This program is in the scripts_figure_5.tar.gz archive. It is   
		  not part of the wrapper and has to be launched sepparately.  
		  The data needed to generate the figures are contained in the   
		  scripts_figure_5/data/ folder. It contains the DAP-seq   
		  peaks of SEP3 and SEP3AG (See ### 2.Peak Calling), the ChIP-seq  
		  peaks of SEP3 and AG (downloaded according to the MATERIAL AND   
		  METHODS section), the Genes regulated by SEP3, AG (downloaded   
		  according to the MATERIAL AND METHODS section) and the genes   
		  either up-regulated or down regulated by both SEP3 and AG   
		  (DEG_SEP3_and_AG.csv). The file all_genes_extended.bed contains  
		  the position of all the genes with an extension of 1500 upstream   
		  the TSS and 600bp downstream the TTS.  
</pre>  
usage : To launch the code, follow these instructions :  
>	$ tar -xzf  scripts_figure_5.tar.gz  
>	$ cd scripts_figure_5  
>	$ cd scripts  
>	$ ./run_all.sh  
   
### 11 Comparing binding sites quality in SEP3AG DAP-seq peaks and SEP3/AG ChIP-seq peaks (figure S8)  
<pre>  
general infos : This program is in the scripts_supp_fig_8.tar.gz archive. It is   
		  not part of the wrapper and has to be launched sepparately.  
		  The data needed to generate the figures are contained in the   
		  scripts_figure_5/data/ folder. It contains the DAP-seq   
		  peaks of SEP3AG (See ### 2.Peak Calling), the ChIP-seq  
		  peaks of SEP3 and AG (downloaded according to the MATERIAL AND   
		  METHODS section), the A.thaliana genome (tair10.fas) and the  
		  SEP3AG TFFM (tffm_first_order.xml).  
</pre>  
usage : To launch the code, follow these instructions :  
>	$ tar -xzf  scripts_supp_fig_8.tar.gz  
>	$ cd scripts_fig_8  
>	$ cd scripts  
>	$ ./run_all.sh  
  
  
  

