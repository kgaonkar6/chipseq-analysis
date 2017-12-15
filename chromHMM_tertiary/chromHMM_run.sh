#!/bin/bash
usage ()
{
cat << EOF
##      reorder - vcf file -w working directory -n case# -t trio[dad,mom,child] -g genelist
##      -i      <INPUT_DIR for bams>
##      -w      <working directory>
##      -n      <state number>
##      -h      <help>
##      -f      <fold change for binarization>
##      -r      <ref name hg19/mm10>
##      -c      < cellmark file>
##	-u 	<User provided coordinates to create overlap enrichement plot; must be a folder with all coordinates eg. /data2/tertiary/Shah_Vijay_vhs01/s205130.UNC_Liver/chromHMM_merged/ChromHMM_model_10/COORDS/hg19/>
EOF
}


#set -x

echo "Options specified: $@"

while getopts ":t:w:n:h:b:r:e:u" OPTION; do
  case $OPTION in
        h) usage
        exit ;;
        i) INPUT_DIR=$OPTARG ;;
        w) WORK_DIR=$OPTARG ;;
        n) state_num=$OPTARG ;;
        f) fold=$OPTARG ;;
        r) ref=$OPTARG ;;
        c) cellmark=$OPTARG ;;
	u) USER_COORDS=$OPTARGS ;;
        \?) echo "Invalid option: -$OPTARG. See output file for usage."
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage."
       usage
       exit ;;
  esac
done




echo "Number of states required" $state_num
echo "Reference used" $ref
echo "Fold change used" $fold
echo "Input directory used" $INPUT_DIR
echo "input file with cell/mark information" $cellmark
echo "Work directory" $WORK_DIR

# Runs standard chromHMM run. Please feel free to try out other paramters for chromHMM. 
java -mx6G -jar /data5/bsi/bictools/src/chromhmm/1.12/ChromHMM.jar BinarizeBam -f $fold /data5/bsi/bictools/src/chromhmm/1.12/CHROMSIZES/${ref}.txt $INPUT_DIR $cellmark $WORK_DIR/Binarized_Bam_fc${fold} 
java -mx6G -jar /data5/bsi/bictools/src/chromhmm/1.12/ChromHMM.jar LearnModel $WORK_DIR/Binarized_Bam_fc${fold} $WORK_DIR/ChromHMM_model_${state_num} ${state_num} $ref

# To obtain overlap enrichment plots for user defined coordinate files. 
# Enhancer bed file from is reguralrly used to annotate enhancer regions. Refer to data2/tertiary/Shah_Vijay_vhs01/s205130.UNC_Liver/chromHMM_merged/ChromHMM_model_10/COORDS/hg19/
# You can also provide RNA data to obtain enrichment plots

if [[ $USER_COORDS != "" ]];then
for file in $(ls $WORK_DIR/ChromHMM_model_${state_num}/*segments.txt);do
        java -mx6G -jar /data5/bsi/bictools/src/chromhmm/1.12/ChromHMM.jar OverlapEnrichment $file $USER_COORDS $(basename $file _segments.txt)
        done
fi
