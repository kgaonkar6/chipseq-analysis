#!/bin/bash

usage ()
{
cat << EOF
##      run_chipTertiary.sh -d <input folder> -c <run_info.txt> 
##	-d	Provide path to secondary folder with "/macs2out" "/sicer" "/mapout" >
##	-c 	Provide full path to run info config file
##	-h	help	
EOF
}


#set -x

echo "Options specified: $@"

while getopts ":d:c:h:" OPTION; do
  case $OPTION in
        h) usage
        exit ;;
	f) input_folder=$OPTARG ;;
	c) runinfo=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG. See output file for usage."
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage."
       usage
       exit ;;
  esac
done

if [[ -z "$Bedfolder" || ! -f $runinfo ]]
then
        echo -e "Must provide at least required options:\n1) Bed files in a folder\n2)run_info.txt\n"
        usage
        exit 1;
fi

source $runinfo

if ($analysis == "diffBind" )
then
source $diffBind_info
for mark in $marks_list
do
#Merge the bed files according to the peak caller with the highest pvalue/qvalue and also get individual bed files for diffBind
qsub -N "MergeBed_${mark}" mergeBed.sh $diffBind_info
done

#Run diffBind for the marks from the config file and assign the regions to genes
qsub -N "Samplesheet_diffBind" -hold_jid "MergeBed_${mark}" run_diffBind.sh $input_dir/ $mark $output_dir/diffBind_samplesheet.txt
fi

#if ( $analysos == "ngsplot")
#then
#ngsplot on given bam and bed file
#qsub -N "Ngsplot" $code/run_NGSpplot.sh $Nfsplot_info
#fi


#if ( $analysis == "ChromHMM")
#then
#Run ChromHMM to identify the functional enrichment of the given marks across the genome
#qsub -N "ChromHMM default" $code/run_ChromHMM.sh $ChromHMM_info
#fi

#if ($analysis == "ROSE")
#then
#Run ROSE to detect super enhancer
#qsub -N "ROSE_default" $code/run_ROSE.sh $ROSE_info
#fi

#if ($analysis == "TFBS")
#then
#run MEME/HOMER for tfbs prediction
#qsub -N "TFBS_prediction" $code/run_TFBS.sh $TFBS_info
#fi



#set +x
