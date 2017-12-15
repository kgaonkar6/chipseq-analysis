#! /usr/bin/env bash

usage ()
{
cat << EOF
##      
##      -i      <list of full paths to sample bam file>
##      -w      <working directory>    
##      -f      <freq>      
##      -o      <output file name>
##      -b      <bed file>
##      -a      <Y/TRUE to run bior annotate>
##	-h	<help>    
EOF
}


#set -x

echo "Options specified: $@"

while getopts ":i:w:f:o:b:a:h:" OPTION; do
  case $OPTION in
        h) usage
        exit ;;
        i) info_file=$OPTARG ;;
        w) output_folder=$OPTARG ;;
        f) freq=$OPTARG ;;
        o) output_file=$OPTARG ;;
	b) bed_file=$OPTARG ;;
        a) run_bior=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG. See output file for usage."
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage."
       usage
       exit ;;
  esac
done

ref="/data2/bsi/reference/sequence/human/ncbi/hg19/allchr.fa"


if [[ ! -f $info_file && $freq == "" ]]
then
echo " Please provide the following parameters to run_varscan.sh"
usage
exit 1;
fi 

VarScan_path="/data2/bsi/staff_analysis/m078940/variant_caller/varscan/varscan-master"
Java_path="/usr/local/biotools/java/jdk1.7.0_03/bin"

list=$(for file in $(cat $info_file);do echo $file;done)
#sample list for vcf header
awk -F "/" '{print $NF}' $info_file > $output_folder/sample_name.txt
sed -i 's/.bam//g' $output_folder/sample_name.txt

#pileup
/data5/bsi/bictools/alignment/samtools/samtools-1.2/samtools mpileup -q 20 --ff DUP -B -f $ref -b $info_file > $output_folder/out.pileup

#varscan
$Java_path/java -jar $VarScan_path/VarScan.v2.4.1.jar mpileup2snp $output_folder/out.pileup --vcf-sample-list $output_folder/sample_name.txt --output-vcf --min-var-freq $freq --min-reads 2 > $output_folder/$output_file
bior_input=$output_folder/$output_file

#subset variants
if [[ -f $bed_file ]]
then
name=$(basename $output_folder/$output_file .vcf)
grep "^#" $output_folder/$output_file > $output_folder/${name}_subset.vcf
bedtools intersect -wa -a $output_folder/$output_file -b $bed_file >> ${name}_subset.vcf
bior_input=$output_folder/${name}_subset.vcf
fi

#run bior
BIOR_PATH="/data5/bsi/bictools/scripts/bior_annotate/v2.4.1/"
BIOR_CATALOGS="/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0/annotation/catalog_file"
BIOR_DRILLS="/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0/annotation/drill_file"
TOOL_INFO="/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0.2/config/tool_info.txt"
MEMORY_INFO="/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0.2/config/memory_info.txt"



if [[ $run_bior == "Y" || $run_bior == "TRUE" ]]
then
$BIOR_PATH/bior_annotate.sh -v $bior_input -O $output_folder -o ${name}_bior_annotate -c $BIOR_CATALOGS -d $BIOR_DRILLS -T $TOOL_INFO -M $MEMORY_INFO -t 1
fi



