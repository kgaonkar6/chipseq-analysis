#!/bin/bash
usage ()
{
cat << EOF
##      reorder -f state ordering file -n state number -w working directory -b path to binarized bam -t tss anchor file -e tes anchorfile -r ref 
##      -t      <TSS anchorfile>
##      -w      <working directory>
##      -n      <state number>
##      -h      <help>
##      -b      <Binarized data>
##      -r      <ref name hg19/mm10>
##	-e	<TES anchorfile>	
##	-f 	< state ordering file>
EOF
}


#set -x

echo "Options specified: $@"

while getopts ":t:w:n:h:b:r:e:f:" OPTION; do
  case $OPTION in
        h) usage
        exit ;;
        t) TSS=$OPTARG ;;
        w) WORK_DIR=$OPTARG ;;
        n) state_num=$OPTARG ;;
        b) BINARIZEDBAM=$OPTARG ;;
        r) ref=$OPTARG ;;
	e) TES=$OPTARG ;;
	f) stateordering=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG. See output file for usage."
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage."
       usage
       exit ;;
	h) echo "Option -$OPTARG requires an argument. See output file for usage."
       usage
       exit ;;
  esac
done

java -mx6G -jar /data2/bsi/staff_analysis/m139467/ChromHMM/edu/mit/compbio/ChromHMM/ChromHMM/ChromHMM.jar Reorder -o $stateordering $WORK_DIR/model_${state_num}.txt $WORK_DIR
java -mx6G -jar /data2/bsi/staff_analysis/m139467/ChromHMM/edu/mit/compbio/ChromHMM/ChromHMM/ChromHMM.jar MakeSegmentation $WORK_DIR/model_${state_num}.txt $BINARIZEDBAM $WORK_DIR
for segment_bed in $(ls $WORK_DIR/*segments.bed);do java -mx6G -jar /data2/bsi/staff_analysis/m139467/ChromHMM/edu/mit/compbio/ChromHMM/ChromHMM/ChromHMM.jar OverlapEnrichment $segment_bed COORDS/$ref $(basename $segment_bed _segments.bed)_enhancer;done
for segment_bed in $(ls $WORK_DIR/*segments.bed);do java -mx6G -jar /data2/bsi/staff_analysis/m139467/ChromHMM/edu/mit/compbio/ChromHMM/ChromHMM/ChromHMM.jar NeighborhoodEnrichment $segment_bed $TSS $(basename $segment_bed _segments.bed)_TSS;done
for segment_bed in $(ls $WORK_DIR/*segments.bed);do java -mx6G -jar /data2/bsi/staff_analysis/m139467/ChromHMM/edu/mit/compbio/ChromHMM/ChromHMM/ChromHMM.jar NeighborhoodEnrichment $segment_bed $TES $(basename $segment_bed _segments.bed)_TES;done
