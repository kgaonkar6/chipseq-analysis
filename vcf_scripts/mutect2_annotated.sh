#!/bin/bash
wd=$1
out=$2

cd $wd
for file in $(ls chr*.vcf.gz)
do
chr=$(basename $file .vcf.gz)
zcat ${chr}.vcf.gz|grep -v "##SAMPLE=<ID" |sed 's/VCFv4.2/VCFv4.1/g' > ${chr}_v4.1.vcf
if [[ -s ${chr}_v4.1.vcf ]] 
then
bgzip ${chr}_v4.1.vcf
tabix -p vcf ${chr}_v4.1.vcf.gz
fi
done

cat_in=$(echo $(for file in $(ls *v4.1.vcf.gz);do echo "-v" $file ;done))
echo $cat_in
/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0.1/CatVariants.sh  $cat_in -t ../../../docs/config/tool_info.txt -m ../../../docs/config/memory_info.txt -o ${out}.vcf
cat ${out}.vcf | perl /data2/bsi/staff_analysis/m078940/projects/braggio/fixDP.pl > ${out}_DP.vcf

bgzip ${out}_DP.vcf
tabix -p vcf ${out}_DP.vcf.gz

/data5/bsi/bictools/scripts/bior_annotate/v2.8.1/trunk/bior_annotate.sh -o ${out}_annot -v ${out}_DP.vcf.gz -O /data2/delivery/Sekulic_Aleksandar_sekulic/161014_SN730_0470_BC9LF0ACXX/secondary/AK_Melanoma/variants/Mutect2/AK_Melanoma/ -T ../../../docs/config/tool_info.txt -M ../../../docs/config/memory_info.txt -t 1 -c /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/4.0.1/annotation/catalog_file.hg38 -L -d /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/4.0.1/annotation/drill_file.hg38


