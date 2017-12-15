#!/bin/bash
#get variant in table form

#$ -N Variant2Table
#$ -m a
#$ -M Gaonkar.Krutika@mayo.edu
#$ -q 1-day
#$ -cwd
#$ -l h_vmem=30G
#$ -l h_stack=10M
#$ -pe threaded 4
#$ -notify 


echo "-i vcf file -c columns -d output directory -t tool_info -f format sections"

while getopts i:c:d:t:f: OPTION; do
        case $OPTION in
                i) vcf=$OPTARG;;
                c) columns_field=$OPTARG;;
		f) columns_format=$OPTARG;;
                d) output_dir=$OPTARG;;
                t) TOOL_INFO=$OPTARG;;
                :)echo "check parameters"
                exit;;

        esac
done
echo "read all the arguments"

#parameters
GATK=$(grep -w 'GATK' $TOOL_INFO| awk -F \= '{print $2}')
echo $GATK
REF=$(grep 'REF_GENOME' $TOOL_INFO| awk -F \= '{print $2}')
echo $REF
JAVA=$(grep 'JAVA7' $TOOL_INFO| awk -F \= '{print $2}')
echo $JAVA

fields=""
  
for col in $(echo $columns_field|tr "," " ") 
do 
		fields=$fields" -F "$col
done

format=""
echo $columns_format
for format in $(echo $columns_format|tr "," " ")
do
	format_fields=$format_fields" -GF "$format
done

name=$(basename $vcf)

echo -e "$JAVA/java -XX:CompileThreshold=1000 -XX:ReservedCodeCacheSize=128m -Xmx10g -Xms5g -Djava.io.tmpdir=$output_dir/temp/ -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable -R $REF -V $vcf $fields $format_fields -o $output_dir/${name}.table"
$JAVA/java -XX:CompileThreshold=1000 -XX:ReservedCodeCacheSize=128m -Xmx10g -Xms5g -Djava.io.tmpdir=$output_dir/temp/ -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable -R $REF -V $vcf $fields $format_fields -o $output_dir/${name}.table











