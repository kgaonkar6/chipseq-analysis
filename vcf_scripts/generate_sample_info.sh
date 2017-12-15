#!/bin/sh
# generate TCGA sample_info with split fastqs

if [ $# != 3 ];
then
	echo "usage: <metadata txt file> <fastq dir> <output sample_info.txt>";
	exit 1;
fi				
#set -x
metadata=$1
fastq_dir=$2
sample_info=$3

cd $fastq_dir
col_array=(ID NB NT TM TP)
samplenames=""
groupnames=""
echo "" > $sample_info

for sample in $(cut -f1 $metadata | sed 1d)
do
	pair=""
	for i in $(seq 2 5);
	do
		bam_column=$(grep $sample $metadata | cut -f$i )
		bam_array=($(echo $bam_column | tr "," " "))
		for j in $(seq ${#bam_array[@]})
		do
			bam=${bam_array[$j-1]}
			bamname=$(basename $bam)
			if [ "$bam" != "/" ]
			then
				fastqs=$(ls $bamname.R1.fastq.gz.* $bamname.R2.fastq.gz.* | sort -t"." -k6 | tr "\n" "\t")
				if [ "$j" == "1" ]
				then
					echo -ne "\n#FASTQ:${sample}_${col_array[$i-1]}=${fastqs}" >> $sample_info
					samplenames="${samplenames}${sample}_${col_array[$i-1]}:"
					if [ "$pair" == "" ]
					then
						pair="${sample}_${col_array[$i-1]}"
					else
						pair="${pair} ${sample}_${col_array[$i-1]}"
					fi
				else
					echo -ne "\t${fastqs}" >> $sample_info
				fi
			fi
		done
	done
	modified_sample=$(echo $sample | tr "-" "_")
	echo -e "\n${modified_sample}=\"${pair}\"" >> $sample_info
	echo "" >> $sample_info
	groupnames="${groupnames}${modified_sample}:"
done


echo -e "SAMPLENAMES:\n${samplenames}\n"
echo -e "GROUPNAMES:\n${groupnames}\n"


