#!/bin/bash
work_dir=$1
mark=$2
file_list=$3
SEQ_DIR=$4

#cd $work_dir
#for file in $(ls ./mapout/*$mark*); do name=$(echo $(basename $file)|cut -d "." -f1|cut -d "H" -f1|sed 's/_$//');input=$(echo ./mapout/$name*input*);echo -e $name,1,$file,${name}_input,$input,$(echo ./peakout/$name*$mark*),bed;done

echo "SampleID,Condition,bamReads,ControlID,bamControl,Peaks,PeakCaller"
for name in $(cat $file_list|cut -d "K" -f1|sed 's/_$//');do bam=$(echo $SEQ_DIR/mapout/$name*$mark*s1.bam);input=$(echo $SEQ_DIR/mapout/$name*input*);bed=$(echo $WORK_DIR/mergeBed/best_$name*$mark*);group=$(grep "$name" $file_list|cut -f2);echo -e $name\,$group\,$bam\,${name}_input\,$input\,$bed\,bed;done

