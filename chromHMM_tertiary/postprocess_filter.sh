#/bin/bash



echo $1

WORK_DIR=$1
segment_file=$(echo $(ls $1/*segments.bed))

names=$(echo $(for file in $(ls $1/*segments.bed);do basename $file;done))
#obtain uninobed for chromhmm segment files
echo -e "chr\tstart\tstop\t"$names|tr " " "\t" > $WORK_DIR/unionbed_segments.txt

bedtools unionbedg -i $(echo $segment_file) >> $WORK_DIR/unionbed_segments.txt

awk -F "\t" '{if($1!="chrY" || $1!="chrX" || $1== "chr") print $0}' $WORK_DIR/unionbed_segments.txt > $WORK_DIR/unionbed_segments_nosexchr.txt










