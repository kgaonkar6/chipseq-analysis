#!/bin/bash

snp_file=$1
GEO_list=$2
WORK_DIR=$3

mkdir $WORK_DIR
mkdir $WORK_DIR/intersectBed_out
mkdir $WORK_DIR/snp_GEO

CODE_DIR=$(dirname $0)
cat $GEO_list|while read LINE;do type=$(echo $LINE|cut -d " " -f3);inputfile=$(echo $LINE|cut -d " " -f2);GEOid=$(echo $LINE|cut -d " " -f1);ref=$(echo $LINE|cut -d " " -f4);
echo $GEOid
if [[ $ref == "hg19" && $type == "ChromHMM" ]]
then
sortBed -i $inputfile|intersectBed -a $snp_file -b stdin -wao > $WORK_DIR/intersectBed_out/$(basename $snp_file .bed)_$(basename $inputfile .bed).bed
for state in CTCF CTCF+Enhancer CTCF+Promoter Enhancer Heterochromatin Poised_Promoter Promoter Repeat Repressed Transcribed 
do
echo -e "SNPid\t$(basename $inputfile .bed)_${state}" > $WORK_DIR/snp_GEO/$(basename $snp_file .bed)_$(basename $inputfile .bed)_ChromHMM_${state}.bed
$CODE_DIR/snptopeakadd.pl $snp_file $WORK_DIR/intersectBed_out/$(basename $snp_file .bed)_$(basename $inputfile .bed).bed $WORK_DIR/snp_GEO/$(basename $snp_file .bed)_$(basename $inputfile .bed)_ChromHMM_${state}.bed $state
done
fi
if [[ $ref == "hg19" && $type == "bed" ]]
then
sortBed -i $inputfile|intersectBed -a $snp_file -b stdin -c > $WORK_DIR/intersectBed_out/$(basename $snp_file .bed)_${GEOid}.bed
echo -e "SNPid\t$(basename $snp_file .bed)_$GEOid" > $WORK_DIR/snp_GEO/$(basename $snp_file .bed)_${GEOid}.bed
$CODE_DIR/snptopeak.pl $snp_file $WORK_DIR/intersectBed_out/$(basename $snp_file .bed)_${GEOid}.bed $WORK_DIR/snp_GEO/$(basename $snp_file .bed)_${GEOid}.bed
fi
done

if  [[ $(ls $WORK_DIR/snp_GEO/*|wc -l) > 1 ]]
then
echo -e "paste <( cut -f1 \$(ls $WORK_DIR/snp_GEO/*|head -1))" $(for file in $(ls $WORK_DIR/snp_GEO/*);do echo "<(cat $file|cut -f2)";done) "> $WORK_DIR/combined_snp_peaksinfo.txt" > $WORK_DIR/merge_snpGEO.job
source $WORK_DIR/merge_snpGEO.job
else
cp $(ls $WORK_DIR/snp_GEO/*|head -1) $WORK_DIR/combined_snp_peaksinfo.txt
fi


