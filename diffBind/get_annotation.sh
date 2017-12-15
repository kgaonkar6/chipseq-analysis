#!/bin/bash
### give diffBind report folder as input
diffBind_report=$1
cd $diffBind_report
for file in $(ls *.txt); do sed 's/ /\t/g' $file|sed 's/"//g'|awk -F "\t" '{if($1!="Chr") print $2"\t"$3"\t"$4"\tPeak"NR"\t"$8}'> $(basename $file .txt).bed ; /data2/bsi/staff_analysis/m139467/HOMER/bin/annotatePeaks.pl $(basename $file .txt).bed hg19 -annStats $(basename $file .txt).annoStats|sed  's/Peak\ Score/log2FC/g'> $(basename $file .txt).anno.xls;done

