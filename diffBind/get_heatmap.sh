#! /bin/bash

#case_num=$1

if [[ $1 == "" || $1 == "-h" ||$2 == ""|| $3 == "" ]]
then
echo -e "\n##\tPlease provide required files as input parameter\t##\n"
echo -e "##\tUsage: ./get_heatmap <metadata file> <DEup_gene.pdf file> <DEdown_gene.pdf file>\t\t##\n"
exit 1;
fi

file=$1
DEup_gene=$2
DEdown_gene=$3


##file="SL1_${1}_6marks.raw.norm.RNAseq.rank.final.txt"

name=$(basename $file .txt)

awk -F '\t' '{if ($21 >= 1 && $22 == "TRUE") print $8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18}' $file | awk '{for (i=1;i<=NF;i++) if ($i >=1) $i=1} {print $0}' > ${name}.up.txt
awk -F '\t' '{if ($21 <= -1 && $22 == "TRUE") print $8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18}' $file | awk '{for (i=1;i<=NF;i++) if ($i >=1) $i=1} {print $0}' > ${name}.down.txt 
 
sort -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n ${name}.up.txt |uniq -c | sort -k1,1rn > ${name}.up.heatmap.txt 
sort -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n ${name}.down.txt|uniq -c | sort -k1,1rn > ${name}.down.heatmap.txt

##will need path to heatmap.R #

Rscript heatmap.R up ${name}.up.heatmap.txt $DEup_gene
Rscript heatmap.R down ${name}.down.heatmap.txt $DEdown_gene

sed -i '1i Count\tH3K4me3_rank\tH3K4me1_rank\tH3K27ac_rank\tH3K36me3_rank\tH3K9me3_rank\tH3K27me3_rank' ${name}.up.heatmap.txt
sed -i '1i Count\tH3K4me3_rank\tH3K4me1_rank\tH3K27ac_rank\tH3K36me3_rank\tH3K9me3_rank\tH3K27me3_rank' ${name}.down.heatmap.txt

