#! /usr/bin/env bash
trap "exit 100;" ERR
####################
## Script Options ##
####################
usage ()
{
cat << EOF
######################################################################
##	script to run RealignertargetCreator on a input file bam 
##	using GATK RealignerTargetCreator module
##	output file created in output folder named Realigner.intervals
##
## Script Options:
##	-o	<output>	-	(REQUIRED)	path to output directory 
##	-i	<inputbam>	-	(REQUIRED)	input bam full path
##	-t	<TOOL_INFO>	-	(REQUIRED)	full path to tool info file
##	-m	<MEMORY_INFO>	-	(REQUIRED)	full path to memory info file
##	-v	<knownvariants>	-	(optional)	known variants (can be specified multiple times)
##	-r	<region>	-	(optional)	region bed file or chromosome number.
##	-c	<chr>		-	(optional)	which chromosome example : chr4
##	-h			- 			Display this usage/help text (No arg)
#############################################################################
EOF
}
echo "Options specified: $@"

while getopts ":o:i:t:m:v:r:c:h" OPTION; do
  case $OPTION in
	v) known_variants="$known_variants $OPTARG" ;;
	r) region="-L "$OPTARG ;;
	c) chr=$OPTARG ;;
	h) usage
	exit ;;
	o) output=$OPTARG ;;
	i) inputbam="$inputbam $OPTARG" ;;
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
 
if [ ! -d "$output" ] || [ -z "$inputbam" ] || [ ! -s "$TOOL_INFO" ] || [ ! -s "$MEMORY_INFO" ];
then
	echo "Must provide at least required options. See output file for usage." >&2
	usage
	exit 1;
else
	if [ ! -s "$TOOL_INFO" ]
	then
		echo "configuration file : $TOOL_INFO is corrupted" >&2
		exit 1;
	else
		source $TOOL_INFO
		source $WORKFLOW_PATH/shared_functions.sh
		check_file $MEMORY_INFO
		check_dir_exist $output
		for i in $inputbam
		do
			check_file $i
		done	
		if [ "$known_variants" ]
		then
			for i in $known_variants
			do
				check_file $i
			done	
		fi
	fi		
fi		
	
email_param="$@"	
START=$(date +%s)
source $MEMORY_INFO

if [ ! -d $output/temp ]
then
	mkdir -p $output/temp
fi
## local variables
check_variable "$TOOL_INFO:REF_GENOME" $REF_GENOME
check_variable "$TOOL_INFO:GATK" $GATK
check_variable "$TOOL_INFO:JAVA7" $JAVA7
check_variable "$TOOL_INFO:THREADS" $THREADS
check_variable "$TOOL_INFO:RealignerTargetCreator_JVM" $RealignerTargetCreator_JVM
if [ $DEBUG_MODE == "YES" ];then set -x; fi
known=""
if [ "$known_variants" ]
then
	for i in $known_variants
	do
		known=$known"--known $i "
	done
fi
bams=""
for i in $inputbam
do
	bams=$bams"-I $i "
done	


gatk_params="-R $REF_GENOME -et NO_ET -K $GATK/Hossain.Asif_mayo.edu.key "
$JAVA7/java $RealignerTargetCreator_JVM -Djava.io.tmpdir=$output/temp/  \
-jar $GATK/GenomeAnalysisTK.jar \
-nt $THREADS \
-T RealignerTargetCreator \
--fix_misencoded_quality_scores \
-o $output/Realigner.intervals $bams $known $RealignerTargetCreator_params $gatk_params $region         
if [ $? -gt 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/Realigner.intervals -m Realignment.sh -s RealignerTargetCreator.sh -p "$email_param" -l $LINENO -M "failed to execute RealignerTargetCreator from GATK"
	exit 100;
fi	
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "RealignerTargetCreator for $inputbam took $DIFF seconds"
