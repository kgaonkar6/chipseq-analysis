#!/bin/bash 

runinfo=$1

source $runinfo

#mergeBed for peak files with q value< 0.01 in macs and FDR <=0.01
mkdir $WORK_DIR/mergeBed
mkdir $WORK_DIR/mergeBed/tmp
mkdir $WORK_DIR/mergeBed/delivery


if [ "$pkcalled_macs2" == "yes" ]
then
caller="macs"
for mark in $marks_macs2
do
echo $mark
for file in $SEC_DIR/macs2out/*$mark*encodePeak
do
perl /data2/bsi/staff_analysis/m139467/diffBind/get_format_bed.pl $file $caller $WORK_DIR/mergeBed
done
cat $SEC_DIR/macs2out/*$mark*encodePeak |sortBed -i stdin| $BEDTOOLS/mergeBed -i stdin -c 9 -o max -d $macs2_merge_dist > $WORK_DIR/mergeBed/mergeBed_${mark}.bed
for file in $WORK_DIR/mergeBed/*$mark*formated.bed
do
perl /data2/bsi/staff_analysis/m139467/diffBind/get_intersect_bestv1.1.pl $file $caller $WORK_DIR/mergeBed/mergeBed_${mark}.bed $WORK_DIR/mergeBed
done
done
fi

if [ $pkcalled_sicer == "yes" ]
then
caller="sicer"
for mark in $marks_sicer
do
echo $mark
for file in $SEC_DIR/sicerout/*$mark*islands-summary-FDR1E*
do
perl /data2/bsi/staff_analysis/m139467/diffBind/get_format_bed.pl $file $caller $WORK_DIR/mergeBed
done
cat $WORK_DIR/mergeBed/*$mark*formated.bed |sortBed -i stdin| $BEDTOOLS/mergeBed -i stdin -c 4 -o max -d $sicer_merge_dist >$WORK_DIR/mergeBed/mergeBed_${mark}.bed
for file in $WORK_DIR/mergeBed/*$mark*formated.bed
do
perl /data2/bsi/staff_analysis/m139467/diffBind/get_intersect_bestv1.1.pl $file $caller $WORK_DIR/mergeBed/mergeBed_${mark}.bed $WORK_DIR/mergeBed
done
done
fi


if [ $get_density == "yes" ]
then
for mark in $marks_macs2 $marks_sicer
do
echo "SampleID,Condition,bamReads,ControlID,bamControl,Peaks,PeakCaller" >$WORK_DIR/mergeBed/diffBind_${mark}_samplesheet.txt


for input in $SEC_DIR/mapout/*[Ii]nput*s1.bam
do
sample=$(echo $(basename $input)|awk '{split($0,a,"-[Ii]nput"); print a[1]}')
total=$(samtools view -c $input)
samtools view $input |awk 'BEGIN {FS="\t"; OFS="\t"} {if ((($2 ==83) || ($2 ==147)) && ($9 >0)) print $3,($4-1),($4-1+$9),".","1","-"; else if ((($2 ==83) || ($2 ==147)) && ($9 <0)) print $3,($4-1),($4-1-$9),".","1","-"; else if ((($2 ==99) || ($2 ==163)) && ($9 >0)) print $3,($4-1),($4-1+$9),".","1","+"; else if ((($2 ==99) || ($2 ==163)) && ($9 <0)) print $3,($4-1),($4-1-$9),".","1","+"}' |$BEDTOOLS/intersectBed -a $WORK_DIR/mergeBed/mergeBed_${mark}.bed -b stdin -c |awk -v Tcount=$total '{$NF=($NF*10000000)/Tcount; print $0}'> $WORK_DIR/mergeBed/tmp/$(echo $(basename $input)|cut -d "." -f1 )_${mark}_density.bed

ipfile=$SEC_DIR/mapout/$sample*$mark*s1.bam
total=$(samtools view -c $ipfile)
samtools view $ipfile |awk 'BEGIN {FS="\t"; OFS="\t"} {if ((($2 ==83) || ($2 ==147)) && ($9 >0)) print $3,($4-1),($4-1+$9),".","1","-"; else if ((($2 ==83) || ($2 ==147)) && ($9 <0)) print $3,($4-1),($4-1-$9),".","1","-"; else if ((($2 ==99) || ($2 ==163)) && ($9 >0)) print $3,($4-1),($4-1+$9),".","1","+"; else if ((($2 ==99) || ($2 ==163)) && ($9 <0)) print $3,($4-1),($4-1-$9),".","1","+"}' |$BEDTOOLS/intersectBed -a $WORK_DIR/mergeBed/mergeBed_${mark}.bed -b stdin -c |awk -v Tcount=$total '{ $NF=($NF*10000000)/Tcount; print $0}'> $WORK_DIR/mergeBed/tmp/$(echo $(basename $ipfile)|cut -d "." -f1 )_density.bed

### Binary density
#####echo -e "chr\tstart\tstop\t${sample}_density">$WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed
#####paste $WORK_DIR/mergeBed/tmp/$sample*$mark*density.bed $WORK_DIR/mergeBed/tmp/$sample*[Ii]nput*$mark*density.bed|awk '{if ($5>=$10) print $1"\t"$2"\t"$3"\t"$5-$10;else print  $1"\t"$2"\t"$3"\t0"}'>> $WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed

# Whole number for density#echo -e "chr\tstart\tstop\t${sample}">$WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed
echo -e "chr\tstart\tstop\t${sample}">$WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed
paste $WORK_DIR/mergeBed/tmp/${sample}*$mark*density.bed $WORK_DIR/mergeBed/tmp/$sample*[Ii]nput_$mark*density.bed|awk '{print $1"\t"$2"\t"$3"\t"$5-$10}'|awk '{$NF=sprintf("%.0f",$NF)}{print $1"\t"$2"\t"$3"\t"$4}'>> $WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed
paste $WORK_DIR/mergeBed/tmp/${sample}*$mark*density.bed $WORK_DIR/mergeBed/tmp/$sample*[Ii]nput_$mark*density.bed|awk '{print $1"\t"$2"\t"$3"\t"$5-$10}'|awk '{$NF=sprintf("%.0f",$NF)}{print $1"\t"$2"\t"$3"\t"$4}'>> $WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed

#if ($get_cluster=- "unsupervised")
#then
#echo $sample condition $ipfile ${sample}_Input $input $WORK_DIR/mergeBed/${sample}_${mark}_IP_INPUT_norm_density.bed bed >> $WORK_DIR/mergeBed/diffBind_${mark}_samplesheet.txt
#else
#$CODE_DIR/get_samplesheet.sh $WORK_DIR $mark $list $SEQ_DIR
#fi

echo -e $sample $mark"\t"$(basename $ipfile)"\t"$(basename $input) >> $WORK_DIR/mergeBed/ChromHMM_${mark}_samplesheet.txt
done

$BEDTOOLS/sortBed -i $WORK_DIR/mergeBed/mergeBed_${mark}.bed > $WORK_DIR/mergeBed/mergeBed_${mark}_sorted.bed
echo "paste $(for file in $WORK_DIR/mergeBed/*$mark*formated.bed;do name=$(basename $file .bed) echo "<($BEDTOOLS/intersectBed -a $WORK_DIR/mergeBed/mergeBed_${mark}_sorted.bed -b $file -c| awk '{if (\$NF >1) \$NF=1;print \$NF}')" ;done)|awk '{for(i=1;i<=NF;i++) t+=\$i; print t; t=0}'|awk 'BEGIN{j=1}{print \"${mark}_peak_\"j\""\\t"\"\$1;j++}'> $WORK_DIR/mergeBed/mergeBed_${mark}_sorted_multi_count.bed"|tr "\n" " "  > $WORK_DIR/mergeBed/multiIntersetBed_${mark}.cmd
source $WORK_DIR/mergeBed/multiIntersetBed_${mark}.cmd
sed -i '1i PeakID\t#samples' $WORK_DIR/mergeBed/mergeBed_${mark}_sorted_multi_count.bed


sed -i '1i chr\tstart\tstop\t-log10(FDR)' $WORK_DIR/mergeBed/mergeBed_${mark}_sorted.bed
echo "paste <( cat $WORK_DIR/mergeBed/mergeBed_${mark}_sorted.bed|cut -f-3) <(cat $WORK_DIR/mergeBed/mergeBed_${mark}_sorted_multi_count.bed|cut -f1) $(echo $(for file in $WORK_DIR/mergeBed/*$mark*IP_INPUT_norm_density.bed;do echo "<(cat $file|cut -f4)";done)) <(cat $WORK_DIR/mergeBed/mergeBed_${mark}_sorted.bed|awk '{if(\$NF>0) \$NF=sprintf(\"%.2f\",\$NF)}{print \$NF}') <(cat $WORK_DIR/mergeBed/mergeBed_${mark}_sorted_multi_count.bed|cut -f2)> $WORK_DIR/mergeBed/final_${mark}_IP_INPUT_norm_density.bed" > $WORK_DIR/mergeBed/final_${mark}_IP_INPUT_norm_density.cmd
source $WORK_DIR/mergeBed/final_${mark}_IP_INPUT_norm_density.cmd

#Rscript ${CODE_DIR}/clustering.R $WORK_DIR/mergeBed/final_${mark}_IP_INPUT_norm_density.bed


done
fi




