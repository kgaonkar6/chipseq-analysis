#!/bin/bash
#Select variants for vcf in multisample ggps data

#$ -N SelectVar
#$ -m a
#$ -M Gaonkar.Krutika@mayo.edu
#$ -q 1-day
#$ -cwd
#$ -l h_vmem=30G
#$ -l h_stack=10M
#$ -pe threaded 4
#$ -notify 


echo "-i input vcf -d output directory -r RUN_info file"

while getopts i:d:r: OPTION; do
        case $OPTION in
		i) input_vcf=$OPTARG;;
		r) RUN_INFO=$OPTARG;;
                d) output_dir=$OPTARG;;
                :)echo "check parameters"
                exit;;

        esac
done
echo "read all the arguments"

#parameters
source $RUN_INFO
source $TOOL_INFO
source $SAMPLE_INFO
echo $GATK
echo $REF_GENOME
echo $JAVA7


total_sample=$(grep 'SAMPLENAMES' $RUN_INFO| awk -F \= '{print $2}'|sed 's/"//g')
total_family=$(grep 'GROUPNAMES' $RUN_INFO| awk -F \= '{print $2}'|sed 's/"//g')



for sample in $(echo $total_sample|tr ':' ' ')
do
$JAVA7/java -XX:CompileThreshold=1000 -XX:ReservedCodeCacheSize=128m -Xmx10g -Xms5g -Djava.io.tmpdir=$output_dir/temp/ -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $REF_GENOME -V $input_vcf -sn $sample -env -o $output_dir/persample/$sample.vcf
done


for fam in $(echo $total_family|tr ':' ' ')
do 
fam_sam=$(grep $fam $SAMPLE_INFO|sed 's/s_/-sn s_/g'|cut -d "=" -f2|sed 's/"//g')
$JAVA7/java -XX:CompileThreshold=1000 -XX:ReservedCodeCacheSize=128m -Xmx10g -Xms5g -Djava.io.tmpdir=$output_dir/temp/ -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $REF_GENOME -V $input_vcf $fam_sam -env -o $output_dir/perfamily/$fam.vcf
done






