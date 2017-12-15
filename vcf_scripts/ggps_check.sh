#! /usr/bin/env bash
## Questions/Recommendations/Concerns: Contact Gaonkar.Krutika@mayo.edu

usage ()
{
cat << EOF
##	ggps_check.sh -o secondary output directory -r run_info > output.txt
##      -o      <secondary output directory>    
##      -r      <run_info>      
##      -h      <help>
##	The script checks for errors,intermediate files and final documentation files created by GGPS.
##	If no error exists then the script exists the folder without printing an error message. 
EOF
}
echo "Options specified: $@"

while getopts ":o:r:h" OPTION; do
  case $OPTION in
        h) usage
        exit ;;
        o) output_dir=$OPTARG ;;
        r) run_info=$OPTARG ;;
	\?) echo "Invalid option: -$OPTARG. See output file for usage."
       usage
       exit ;;
    :) echo "Option -$OPTARG requires an argument. See output file for usage."
       usage
       exit ;;
  esac
done

if [[ -z "$output_dir" || -z "$run_info" ]]
then
        echo "Must provide at least required options. See output file for usage."
        usage
        exit 1;
fi

source $run_info
source $TOOL_INFO

#### check for errors in log folders
## if sample number is less than 2 then dont check CNV folders

echo -e "####\tThe script checks for errors,intermediate files and final documentation files created by GGPS.\t####\n####\tIf no error exists then the script exits the folder without printing an error message.\t\t####\n"

echo "Checking for core files in log folder" 
core_files=0
#core files
if [ -f $output_dir/.logs/core.* ]
then
let core_files=$core_files+1
fi
echo "$core_files core files exist."

sample_count=$(echo $SAMPLENAME|tr ":" "\n")

echo "Checking for error files in log folder"
#errors
error_files=0
for error in Error Missing Fatal 
do
grep -i -l -m 1 $error $output_dir/.logs/*.e*
if [[ $DATATYPE != "whole_genome"  && $sample_count>2 ]]
then 
grep -i -l -m 1 $error $output_dir/.tmp/cnv/logs/*.e*
error_files=$(grep -i -l -m 1 $error $output_dir/.tmp/cnv/logs/*.e*|wc -l)
fi
done
echo "$error_files error files exist."
echo "Exiting log folder"


###check for html and readme
echo "Checking for index.html,README,ngsqc.pdf"
README="$output_dir/README"
INDEX="$output_dir/index.html"
NGSQC_PDF="$output_dir/qc/ngsqc/ngsqc.pdf"

for i in $README $INDEX $NGSQC_PDF 
do
if [ ! -f $i ]
then
echo "WARNING:$i missing"
else
echo "Final documentation files are available"
fi

done
 
#####check for vcf filter/ annoatation

echo "Checking for errors and intermediate files in variants folder"
inputvcf="variants.vcf"
outSNVvcf=$inputvcf.SNV.filter.vcf.gz
recalSNVfile=$inputvcf.SNV.vcf.recal.vcf
tranchesSNV=$inputvcf.SNV.vcf.tranches
rscriptSNVFile=$inputvcf.SNV.vcf.Rscript.R
outINDELvcf=$inputvcf.INDEL.filter.vcf.gz
recalINDELfile=$inputvcf.INDEL.vcf.recal.vcf
tranchesINDEL=$inputvcf.INDEL.vcf.tranches
rscriptINDELFile=$inputvcf.INDEL.vcf.Rscript.R
VCF=$output_dir/variants/
VQSR_SNP=$(zcat $VCF/variants.vcf.gz|head -30|grep "VQSRTrancheSNP" |wc -l)
VQSR_INDEL=$(zcat $VCF/variants.vcf.gz|head -30|grep "VQSRTrancheINDEL" |wc -l)
if [[ -f $VCF/$outSNVvcf ||  -f $VCF/$recalSNVfile ||  -f $VCF/$tranchesSNV ||  -f $VCF/$rscriptSNVFile  ||  -f $VCF/$outINDELvcf || -f $VCF/$recalINDELfile ||  -f $VCF/$tranchesINDEL  ||  -f $VCF/$rscriptINDELFile || "$VQSR_SNP" -eq 0 ||"$VQSR_INDEL" -eq 0 ]]
then
echo "The SNV and INDEL filter and merge process was incomplete .Re-run filtervcf.sh"
else 
echo "The SNV and INDEL filter and merge process complete"
fi
chr_count=$(tabix -l $VCF/variants.vcf.gz|wc -l)
chr_runinfo=$(echo $CHRINDEX|tr ":" "\n"|wc -l) 
if [[ "$chr_count" -ne "$chr_runinfo" ]]
then
echo "Missing chromosome check with tabix -l variants.vcf.gz"
fi
echo "Exiting variant folder" 


#####check for NGSQC files
echo "Checking istats and sexck folder"
for group in $(echo $GROUPNAMES|tr ":" "\n")
do
if [[ $group != "NA" ]]
then
#NGSQC_istats_check
FILE1SNV="$output_dir/qc/ngsqc/istats/$group.snv.pseq.gz"
FILE1INDEL="$output_dir/qc/ngsqc/istats/$group.snv.pseq.gz"
FILE2SNV="$output_dir/qc/ngsqc/istats/$group.snv.bychrom.pseq.gz"
FILE2INDEL="$output_dir/qc/ngsqc/istats/$group.snv.bychrom.pseq.gz"
if [[ ! -f $FILE1SNV || ! -f $FILE2SNV || ! -f $FILE1INDEL || ! -f $FILE2INDEL ]]
then
echo "WARNING:istats for $group incomplete. Re-run ngsqc.sh "
fi



#NGSQC_sexck check
FILE1="$output_dir/qc/ngsqc/sexck/$group.sexcheck"
FILE2="$output_dir/qc/ngsqc/sexck/$group.chry.imiss"
FILE3="$output_dir/qc/ngsqc/sexck/$group.chry.lmiss"

if [[ ! -f $FILE1 || ! -f $FILE2 || ! -f $FILE3 ]]
then
echo "WARNING:sexck for $group incomplete. Re-run ngsqc.sh"
fi
fi

if [[ $group == "NA" ]]
then
#NGSQC_istats_check
FILE1SNV="$output_dir/qc/ngsqc/istats/$PROJECTNAME.snv.pseq.gz"
FILE1INDEL="$output_dir/qc/ngsqc/istats/$PROJECTNAME.snv.pseq.gz"
FILE2SNV="$output_dir/qc/ngsqc/istats/$PROJECTNAME.snv.bychrom.pseq.gz"
FILE2INDEL="$output_dir/qc/ngsqc/istats/$PROJECTNAME.snv.bychrom.pseq.gz"
if [[ ! -f $FILE1SNV || ! -f $FILE2SNV || ! -f $FILE1INDEL || ! -f $FILE2INDEL ]]
then
echo "WARNING:istats for $PROJECTNAME incomplete. Re-run ngsqc.sh "
fi


#NGSQC_sexck check
FILE1="$output_dir/qc/ngsqc/sexck/$PROJECTNAME.sexcheck"
FILE2="$output_dir/qc/ngsqc/sexck/$PROJECTNAME.chry.imiss"
FILE3="$output_dir/qc/ngsqc/sexck/$PROJECTNAME.chry.lmiss"

if [[ ! -f $FILE1 || ! -f $FILE2 || ! -f $FILE3 ]]
then
echo "WARNING:sexck for $PROJECTNAME incomplete. Re-run ngsqc.sh"
fi

fi

done
echo "Exiting NGSQC folder"

######check for CNV plots
if [[ $DATATYPE != "whole_genome"  && $sample_count >2 ]]
then
echo "Checking for CNV plots"
CNVwigs=$(ls $output_dir/.tmp/cnv/wigs/)
CNVplots=$(echo $output_dir/qc/cnv/)
for sample in $(echo $SAMPLENAMES|tr ":" "\n")
do
if [[ $(ls $CNVplots/$sample*.png|wc -l) -eq 0 ]]
then
echo "Missing CNV plots for $CNVplots/$sample"
fi
done
echo "Exiting CNV plots"
fi











