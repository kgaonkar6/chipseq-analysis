#!/bin/sh
# split a VCF by sample

if [ $# != 2 ];
then
	echo "usage: <input vcf> <output_dir>";
	exit 1;
fi
	set -x
	input_vcf=$1
	output_dir=$2

for sample in $(zcat $input_vcf | grep "^#CHROM" | cut -f10-)
do

	/usr/local/biotools/java/jdk1.7.0_67/bin/java -jar /data5/bsi/bictools/alignment/gatk/3.3-0/GenomeAnalysisTK.jar -T SelectVariants -R /data2/bsi/reference/sequence/human/ncbi/hg19/allchr.fa -V $input_vcf -sn $sample -o $output_dir/variants.$sample.vcf --excludeNonVariants

	bgzip $output_dir/variants.$sample.vcf
	tabix -p vcf $output_dir/variants.$sample.vcf.gz

done

echo $(date)

