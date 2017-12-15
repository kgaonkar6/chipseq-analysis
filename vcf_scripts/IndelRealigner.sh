#! /usr/bin/env bash

####################
## Script Options ##
####################
usage ()
{
cat << EOF
######################################################################
##	script to run Indelrealigner walker after RealignTragetCreator 
##	walker from GATK on a input file bam
##	output file : realigned.bam in output folder
##
## Script Options:
##	-o	<output>	-	(REQUIRED)	path to output directory 
##	-i	<inputbam>	-	(REQUIRED)	input bam full path
##	-g	<input_intervals>	-	(REQUIRED)	bed file from RealignTargetCreator module
##	-t	<TOOL_INFO>	-	(REQUIRED)	full path to tool info file
##	-m	<MEMORY_INFO>	-	(REQUIRED)	full path to memory info file
##	-v <known_variants>	-	(optional)	variants vcf file; can be specified mltiple times
##	-r	<region>	-		(optional)	region bed file or chromosome number.
##	-h	- Display this usage/help text (No arg)
#############################################################################
EOF
}
echo "Options specified: $@"

while getopts ":o:i:g:t:m:v:r:h" OPTION; do
  case $OPTION in
	v) known_variants="$known_variants $OPTARG" ;;
	r) region="-L "$OPTARG ;;
	h) usage
	exit ;;
	o) output=$OPTARG ;;
	i) inputbam="$inputbam $OPTARG" ;;
    g) input_intervals=$OPTARG ;;
	t) TOOL_INFO=$OPTARG ;;
	m) MEMORY_INFO=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG. See output file for usage." >&2
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage." >&2
       usage
       exit ;;
  esac
done
 
if [ ! -d "$output" ] || [ -z "$inputbam" ] || [ ! -s "$input_intervals" ] || [ ! -s "$TOOL_INFO" ] || [ ! -s "$MEMORY_INFO" ];
then
	echo "Must provide at least required options. See output file for usage." >&2
	usage
	exit 1;
else
	for i in $inputbam
	do
		if [ ! -s $i ]
		then
			echo "BAM file provided : $i doesn't exist" >&2
			usage
			exit 1;
		fi
	done
	if [ "$known_variants" ]
	then
		for i in $known_variants
		do
			if [ ! -s $i ]
			then
				echo "known variant file provided : $i doesn't exist" >&2
				usage
				exit 1;
			fi
		done
	fi	
fi
email_param="$@"	
START=$(date +%s)
source $TOOL_INFO
source $WORKFLOW_PATH/shared_functions.sh
source $MEMORY_INFO
if [ ! -d $output/temp ]
then
	mkdir -p $output/temp
fi

known=""
if [ "$known_variants" ]
then
	for i in $known_variants
	do
		known=$known"--knownAlleles $i "
	done
fi

bams=""
for i in $inputbam
do
	bams=$bams"-I $i "
done	

## local variables
check_variable "$TOOL_INFO:REF_GENOME" $REF_GENOME
check_variable "$TOOL_INFO:GATK" $GATK
check_variable "$TOOL_INFO:JAVA7" $JAVA7
check_variable "$TOOL_INFO:THREADS" $THREADS
check_variable "$TOOL_INFO:IndelRealigner_JVM" $IndelRealigner_JVM
if [ $DEBUG_MODE == "YES" ];then set -x; fi

gatk_params="-R $REF_GENOME -et NO_ET -K $GATK/Hossain.Asif_mayo.edu.key "

$JAVA7/java $IndelRealigner_JVM -Djava.io.tmpdir=$output/temp/ \
-jar $GATK/GenomeAnalysisTK.jar \
-T IndelRealigner \
--out $output/realigned.bam  \
--fix_misencoded_quality_scores \
-targetIntervals $input_intervals $IndelRealigner_params $known $bams $gatk_params $region 
if [ $? -gt 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/realigned.bam -m Realignment.sh -s IndelRealigner.sh -p "$email_param" -l $LINENO -M "failed to execute IndelRealigner from GATK "
	exit 100;
fi
mv $output/realigned.bai $output/realigned.bam.bai  
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "IndelRealigner for $input_bam took $DIFF seconds"
