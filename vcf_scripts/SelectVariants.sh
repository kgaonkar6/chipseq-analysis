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
##	-o	<output>	-	(REQUIRED)	output vcf file full path
##	-i	<input>		-	(REQUIRED)	input vcf full path
##	-t	<TOOL_INFO>	-	(REQUIRED)	full path to tool info file
##	-m	<MEMORY_INFO>	-	(REQUIRED)	full path to memory info file
##	-s	<SN_LIST>	-	(REQUIRED)	-sn <sample1> -sn <sample2> ... 	
##	-h			- 			Display this usage/help text (No arg)
#############################################################################
EOF
}
echo "Options specified: $@"

while getopts ":o:i:t:m:h:s:" OPTION; do
  case $OPTION in
	i) input="$OPTARG" ;;
	h) usage
	exit ;;
	o) output=$OPTARG ;;
	t) TOOL_INFO=$OPTARG ;;
	m) MEMORY_INFO=$OPTARG ;;
	s) SN_LIST=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG. See output file for usage." >&2
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage." >&2
       usage
       exit ;;
  esac
done
 
if [ ! -s "$input" ] || [ ! -s "$TOOL_INFO" ] || [ ! -s "$MEMORY_INFO" ];
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
	fi
fi		
	
START=$(date +%s)
source $MEMORY_INFO

if [ ! -d $PWD/temp ]
then
	mkdir -p $PWD/temp
fi
## local variables
 
echo $SN_LIST;
THREADS=4

gatk_params="-R $REF_GENOME -et NO_ET -K $GATK/Hossain.Asif_mayo.edu.key -env"
$JAVA7/java $RealignerTargetCreator_JVM -Djava.io.tmpdir=$PWD/temp/  \
-jar $GATK/GenomeAnalysisTK.jar \
-nt $THREADS \
-T SelectVariants \
-V $input \
-o $output \
$SN_LIST \
$gatk_params 

#-select "DP > 10"


END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SelectVariants for $inputbam took $DIFF seconds"
