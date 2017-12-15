#!/bin/sh

# ngsplot for TSS +/- xx kb, tested 03/11/2014
# final version, 05/27/2014

###########    parameters need to changed is labeled as <===

###### variable
#SOURCE_DIR=/data5/bsi/bioinf_ext1/s115547.Chip-Seq/references/annotation/hg19
SOURCE_DIR=/data2/bsi/staff_analysis/m139467/
NGSPLOT_PATH=/data5/bsi/bictools/src/ngsplot-2.02/bin   #this is old version  default newest version command this line
#GENOME_VERSION=hg19
GENOME_VERSION=mm10
EXTENSION=5000
EXT=$(awk -v e=$EXTENSION 'BEGIN {print e/1000}')

PLOT_LEVEL=gene
#PLOT_LEVEL=transcript

WORK_DIR=/data2/delivery/Westendorf_Jennifer_mrk4884/150716_SN725_0502_BC77AGACXX/tertiary/ngsplot/5KAroundTSS_down_Cre           #<=====
LOG_FILE=${WORK_DIR}"/"ngsplot.$GENOME_VERSION.$PLOT_LEVEL.$EXTENSION.TSS.log.txt


###### mouse data
#OUT_DIR=/Users/m105265/PI_support_project/Dr_Tamas_Ordog/ChiPseq_130705_SN730_0284_BC29K1ACXX/U12_ngsplot
#BAM_DIR=/Users/m105265/PI_support_project/Dr_Tamas_Ordog/ChiPseq_130705_SN730_0284_BC29K1ACXX/U12_s1_bam
#ALIGNMENT_FILE=${WORK_DIR}"/"temp_mm_bam_list.txt


###### human data
OUT_DIR=/data2/delivery/Westendorf_Jennifer_mrk4884/150716_SN725_0502_BC77AGACXX/tertiary/ngsplot/5KAroundTSS_down_Cre   # <======

## a lit of bam files
BAM_DIR=/data2/delivery/Westendorf_Jennifer_mrk4884/150716_SN725_0502_BC77AGACXX/secondary/chipseq_macs2/mapout       #  <=====
ALIGNMENT_FILE=${WORK_DIR}"/"bam_list.txt

BAM_ARRAY=$( cat $ALIGNMENT_FILE |cut -f 1 |awk -v NAME=${BAM_DIR} '{print NAME"/"$1}' |tr -s "\n" " ")
BAM_FILE=( $( echo $BAM_ARRAY ) )


if [[ ! -d $OUT_DIR ]]; then mkdir $OUT_DIR; fi


if [[ $GENOME_VERSION = "hg19" ]]
then

     if [[ $PLOT_LEVEL = "gene" ]]
      then
      ANNO_FILE=${SOURCE_DIR}"/"refGene.Genes2.hg19.txt

      elif [[ $PLOT_LEVEL = "transcript" ]]
      then
      ANNO_FILE=${SOURCE_DIR}"/"refGene.Transcritps.hg19.txt
      
      fi


elif [[ $GENOME_VERSION = "mm10" ]]
then

     if [[ $PLOT_LEVEL = "gene" ]]
     then
     ANNO_FILE=${SOURCE_DIR}"/"refGene.down.Cre.txt

     elif [[ $PLOT_LEVEL = "transcript" ]]
     then
     ANNO_FILE=${SOURCE_DIR}"/"refGene.Transcritps.mm10.txt

     fi

fi


####### generate TSS file

echo -e "\nStart to generate TSS file from $(basename $ANNO_FILE)" >> $LOG_FILE

FNAME=$(basename ${ANNO_FILE} .txt)

#TAB=`echo -e "\t"`
awk 'BEGIN {FS=OFS="\t"} {if (($1 ~ /^chr[1-9]$/) || ($1 ~ /^chr[1-2][0-9]$/) || ($1 ~ /^chr[M-Y]$/)) print $0}' ${ANNO_FILE} | \
awk 'BEGIN {FS=OFS="\t"} {if ($6 ~ /\+/) print $1,$2,($2+1),$4,$1"/"$2"/"$3"/"$4"/"$5,$6; else if ($6 ~ /\-/) print $1,$3,($3+1),$4,$1"/"$2"/"$3"/"$4"/"$5,$6}' | \
sort -k1,1 -k2,2n -k3,3n > \
${OUT_DIR}"/"$FNAME.TSS.$EXTENSION.bed


####### ngsplot

for ((j=0;j<${#BAM_FILE[@]};j=j+1))
#for ((j=0;j<1;j=j+1))

do

FNAME2=$( basename ${BAM_FILE[j]} .PE.U12.dedup.s1.bam |tr -s "\." "\t" |cut -f 1 )

echo -e "\nStart to generate ngsplot for library $FNAME2" >> $LOG_FILE


## plot   could use sge here ############
$NGSPLOT_PATH"/"ngs.plot.r \
-G $GENOME_VERSION \
-R ${OUT_DIR}"/"$FNAME.TSS.$EXTENSION.bed \
-C ${BAM_FILE[j]} \
-O ${OUT_DIR}"/"${FNAME}.TSS.${EXT}kb_vs_${FNAME2} \
-T ${FNAME}.TSS_${EXT}kb_vs_${FNAME2} \
-F protein_coding -D refseq -I 0 -L $EXTENSION -P 18 \
-GO total -AL bin -FL 150 -MQ 0 -RB 0 -RZ 0 -SC 0,10 -FC 0.02 -FI 0 -H 0


## process output files
cd ${OUT_DIR}
unzip ${OUT_DIR}"/"${FNAME}.TSS.${EXT}kb_vs_${FNAME2}.zip

awk '{if ($0 ~ /me/) print (NR-1)"\t"$1; else printf "%s\t%.4f\n", (NR-1),$1}' \
${OUT_DIR}/${FNAME}.TSS.${EXT}kb_vs_${FNAME2}"/"avgprof.txt | \
head -102 | \
awk 'BEGIN {FS=OFS="\t"} {if (NR==1) print "#Bin","AVE_tag_density"; else print $0}' > \
${OUT_DIR}"/"${FNAME}.TSS.${EXT}kb_vs_${FNAME2}.avgprof.txt


cut -f 3-105 ${OUT_DIR}/${FNAME}.TSS.${EXT}kb_vs_${FNAME2}"/"hm1.txt | \
awk 'BEGIN {FS=OFS="]\t"} {if (NR==1) gsub(/tid/, "#Chr\tstart\tend\tgene_symbol\taccession",$1); if (NR >=2) gsub(/\//, "\t",$1); print}' > \
${OUT_DIR}"/"${FNAME}.TSS.${EXT}kb_vs_${FNAME2}.heatmap.txt


# rm -rf ${OUT_DIR}"/"${FNAME}.mm10.TSS.5kb_vs_${FNAME2}
# rm -rf ${OUT_DIR}"/"${FNAME}.mm10.TSS.5kb_vs_${FNAME2}.heatmap.pdf
# rm -rf ${OUT_DIR}"/"${FNAME}.mm10.TSS.5kb_vs_${FNAME2}.zip

cd $WORK_DIR

echo -e "\nFinish generate ngsplot for library $FNAME2, $(date)\n" >> $LOG_FILE

done

