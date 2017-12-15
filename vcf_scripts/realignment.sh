#! /usr/bin/env bash
trap "exit 100;" ERR

####################
## Script Options ##
####################

usage ()
{
cat << EOF
######################################################################
##	script to run realignemnt on splitted or non splitted bam file with one or more files
##
## Script Options:
##	-d	<input>	-	(REQUIRED)	path to input directory if multiple then sep by ':'
##	-i	<inputbam>	-	(REQUIRED)	full path to input bams, can specify multiple bams sep by ":"
##	-b	<outputbam>	-	(REQUIRED)	name of output BAM file
##	-s	<sample>	-	(REQUIRED)	sample name
##	-o	<output>	-	(REQUIRED)	full path of output folder
##	-r	<run_info>	-	(REQUIRED)	full/path/to/run info file
##	-t	<datatype>	-	(REQUIRED)	exome/whole_genome
##	-c	<chrsplitted>	-	flag to say if BAM is splitted in chromosome; if specified then yes	
##	-f	<recalibration>	-	flag to say if the process involved recalibration. if specified then no
##	-z	<SGE_TASK_ID>	-	SGE_TASK_ID (optional)
##	-h	- Display this usage/help text (No arg)
#############################################################################
EOF
}
echo "Options specified: $@"

while getopts ":d:i:b:s:o:r:t:z:g:cfh" OPTION; do
  case $OPTION in
    g) group=$OPTARG ;;
    z) SGE_TASK_ID=$OPTARG ;;
	c) chrsplitted="yes" ;;
	f) recalibration="no";;
	h) usage
	exit ;;
	d) input=$OPTARG ;;
	i) inputbam=$OPTARG ;;
    b) outputbam=$OPTARG ;;
	s) sample=$OPTARG ;;
	o) output=$OPTARG ;;
	r) run_info=$OPTARG ;;
	t) datatype=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG. See output file for usage." >&2
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage." >&2
       usage
       exit ;;
  esac
done
 
if [ -z "$inputbam" ] || [ -z "$sample" ] || [ -z "$output" ] || [ -z "$run_info" ]  || [ -z "$datatype" ];
then
	echo "Must provide at least required options. See output file for usage." >&2
	usage
	exit 1;
else
	if [ ! -s "$run_info" ]
	then
		echo "configuration file : $run_info is corrupted" >&2
		exit 1;
	else
		source $run_info
		source $TOOL_INFO
		source $WORKFLOW_PATH/shared_functions.sh
		check_file_nonexist $output/$outputbam 
		check_dir_exist $output
		if [[ $datatype != "whole_genome" && $datatype != "exome" ]]
		then
			echo "DATATYPE : $datatype is not defined correctly, it should be whole_genome or exome"
			exit 1;
		fi
	fi	
fi
START=$(date +%s)
email_param="$@"
### if job exits with errortrap
#trap "/bin/rm -rf $output/$outputbam; exit" 1 2 3 6 7 8 9 10 11 12 15 16 17 EXIT

### creating the local variables
check_variable "$run_info:MEMORY_INFO" $MEMORY_INFO
check_variable "$TOOL_INFO:DBSNP_VCF" $DBSNP_VCF
check_variable "$TOOL_INFO:TARGETTED" $TARGETTED
check_variable "$TOOL_INFO:KGENOME_REF" $KGENOME_REF
check_variable "$TOOL_INFO:SAMTOOLS" $SAMTOOLS
check_variable "$TOOL_INFO:DEBUG_MODE" $DEBUG_MODE
check_variable "$TOOL_INFO:PERL" $PERL
if [ $DEBUG_MODE == "YES" ];then set -x; fi

if [ "$chrsplitted" ]
then
	if [ "$SGE_TASK_ID" ]
	then
		chr=`echo $CHRINDEX |tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1`
	fi
	output=$output/chr$chr
	mkdir -p $output
	if [ -f $output/$outputbam ]
	then
		echo "BAM file : $output/$outputbam already exist" >&2 
		usage
		exit 1;
	fi
fi

if [ $datatype == "exome" ]
then
	check_variable "$TOOL_INFO:ONTARGET" $ONTARGET
fi

known_variants=""
### checking for the reference files
if [[ ${#DBSNP_VCF} -ne 0 || $DBSNP_VCF != "NA" ]]
then
	if [ "$chrsplitted" ]
	then
		$TABIX/tabix -h $DBSNP_VCF chr$chr | $TABIX/bgzip -f > $output/dbSNP.vcf.gz
		$TABIX/tabix -f -p vcf $output/dbSNP.vcf.gz
		DBSNP_VCF="$output/dbSNP.vcf.gz"
	fi
	known_variants="-v $DBSNP_VCF" 
fi

if [[ ${#KGENOME_REF} -ne 0 || $KGENOME_REF != "NA" ]]
then
	if [ "$chrsplitted" ]
	then
		$TABIX/tabix -h $KGENOME_REF chr$chr | $TABIX/bgzip -f > $output/KGENOME_REF.vcf.gz
		$TABIX/tabix -f -p vcf $output/KGENOME_REF.vcf.gz
		KGENOME_REF="$output/KGENOME_REF.vcf.gz"
	fi
	known_variants=$known_variants" -v $KGENOME_REF" 
fi


### traversing region
if [ ! -d $output/temp ]
then
	mkdir -p $output/temp
fi

if [[ "$chrsplitted" ]]
then
	if [ $datatype == "whole_genome" ]
	then
		region="chr${chr}"
	else
		if [ $TARGETTED == "YES" ]
		then
			cat $ONTARGET | grep -w chr$chr > $output/target.bed
			if [ `cat $output/target.bed | wc -l` -gt 0 ]
			then
				region="$output/target.bed"
			else
				region="chr${chr}"
			fi	
		else
			region="chr$chr"
		fi	
	fi
else
	if [ $datatype == "exome" ]
	then
		if [ $TARGETTED == "YES" ]
		then
			region="$ONTARGET "
		fi
	fi	
fi	

### make sure all the bam file are present and in GATK acceptable format
j=1
input_bam=""
for i in `echo $input | tr ":" " "`
do
	if [ "$chrsplitted" ]
	then
		i=$i/chr$chr
	fi	
	bam=`echo $inputbam | cut -f "$j" -d ':'`
	if [ ! -s  $i/$bam ]
	then
		$WORKFLOW_PATH/email.sh -f $i/$bam -m realignment.sh -p "$email_param" -l $LINENO
		exit 100;
	fi	
	validate_bam $SAMTOOLS $i/$bam $output
	if [ $? -gt 0 ]
	then
		$WORKFLOW_PATH/email.sh -f $i/$bam -m realignment.sh -p "$email_param" -l $LINENO
		exit 100;
	fi	
	if [ ! -s $i/$bam.bai ]
	then
		$WORKFLOW_PATH/indexbam.sh -i $i/$bam -t $TOOL_INFO
		if [ $? -gt 0 ]
		then
			$WORKFLOW_PATH/email.sh  -f $i/$bam.bai -m realignment.sh -s indexbam.sh -p "-i $i/$bam -t $TOOL_INFO" -l $LINENO
			exit 100;
		fi	
	fi	 
	###check BAM sorting
	sortflag=`$PERL $WORKFLOW_PATH/checkBAMsorted.pl -i $i/$bam -s $SAMTOOLS`
	if [ $sortflag == 0 ]
	then
		if [ -f $output/$bam ]; then rm $output/$bam;fi
		$WORKFLOW_PATH/sortbam.sh -i $i/$bam -o $output/$bam -d $output -r coordinate -t $TOOL_INFO -m $MEMORY_INFO
		if [ $? -gt 0 ]
		then
			$WORKFLOW_PATH/email.sh -f $output/$bam -m realignment.sh -s sortbam.sh -p "-i $i/$bam -o $output/$bam -d $output -r coordinate -t $TOOL_INFO -m $MEMORY_INFO" -l $LINENO
			exit 100;
		fi	
		validate_bam $SAMTOOLS $output/$bam $output 			
	else
		if [ -f $output/$bam ]; then rm $output/$bam;fi
		ln -s $i/$bam $output/$bam
		if [ -f $output/$bam.bai ]; then rm $output/$bam.bai;fi
		ln -s $i/$bam.bai $output/$bam.bai 	
	fi											
	let j=j+1
	input_bam=$input_bam"-i $output/$bam "
done 	

#### 2 modules of GATK to run the process
if [[ "$chrsplitted"  ]]
then
	s_param="-o $output $input_bam $known_variants -t $TOOL_INFO -m $MEMORY_INFO -r $region -c chr$chr" 
else
	if [ "$region" ]
	then
		s_param="-o $output $input_bam $known_variants -t $TOOL_INFO -m $MEMORY_INFO -r $region "
	else
		s_param="-o $output $input_bam $known_variants -t $TOOL_INFO -m $MEMORY_INFO"
	fi
fi	

echo $s_param | xargs /data2/bsi/staff_analysis/m139467/RealignerTargetCreator.sh
if [ $? -gt 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/Realigner.intervals -m realignment.sh -s RealignerTargetCreator.sh -p "$s_param" -l $LINENO
	exit 100;
fi
	
if [[ "$chrsplitted" ]]
then
	s_param="-o $output $input_bam -g $output/Realigner.intervals $known_variants -t $TOOL_INFO -m $MEMORY_INFO -r chr$chr " 
else
	s_param="-o $output $input_bam -g $output/Realigner.intervals $known_variants -t $TOOL_INFO -m $MEMORY_INFO"
fi

echo $s_param | xargs /data2/bsi/staff_analysis/m139467/IndelRealigner.sh
if [ $? -gt 0 ]
then
	$WORKFLOW_PATH/email.sh -f $output/realigned.bam -m realignment.sh -s IndelRealigner.sh -p "$s_param" -l $LINENO
	exit 100;
fi	
validate_bam $SAMTOOLS $output/realigned.bam $output
if [ $? -eq 0 ]
then
	rm $output/Realigner.intervals
	if [ -f $output/target.bed ];then rm $output/target.bed;fi
	for i in $input_bam
	do
		if [ -s $i ];then rm $i;fi
		if [ -s $i.bai ];then rm $i.bai;fi
	done
else
	$WORKFLOW_PATH/email.sh -f $output/realigned.bam -m realignment.sh -s IndelRealigner.sh -p "$s_param" -l $LINENO
	exit 100;	
fi

if [ "$recalibration" ]
then
	$WORKFLOW_PATH/flagstat.sh -m $MEMORY_INFO -i $output/realigned.bam -o $output/flagstat.txt -t $TOOL_INFO -a SAMTOOLS
fi	

### create the idxstats file for BAM file
if [ $DEBUG_MODE == "NO" ]
then
	if [[ "$chrsplitted"  ]]
	then
		if [ $group ]
        then
            s_param="STORE $run_info $output/realigned.bam $group.chr$chr.realign"

        else
            s_param="STORE $run_info $output/realigned.bam $sample.chr$chr.realign"
        fi
    else
        if [ $group ]
        then
            s_param="STORE $run_info $output/realigned.bam $group.realign"
        else
            s_param="STORE $run_info $output/realigned.bam $sample.realign"
        fi
    fi
	echo $s_param | xargs $WORKFLOW_PATH/validate_idxstat.sh
fi	

#### cleanup
if [ -s $output/KGENOME_REF.vcf.gz ]; then rm $output/KGENOME_REF.vcf.gz; fi
if [ -s $output/KGENOME_REF.vcf.gz.tbi ]; then rm $output/KGENOME_REF.vcf.gz.tbi; fi

if [ -s $output/dbSNP.vcf.gz ]; then rm $output/dbSNP.vcf.gz; fi
if [ -s $output/dbSNP.vcf.gz.tbi ]; then rm $output/dbSNP.vcf.gz.tbi; fi


END=$(date +%s)
DIFF=$(( $END - $START ))
echo "realignment for $sample took $DIFF seconds" 
