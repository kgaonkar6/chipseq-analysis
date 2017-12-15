#!/usr/bin/env bash
#exec 3>&2
#exec 2> /dev/null
source /usr/local/biotools/bior_scripts/2.4.0/PKG_PROFILE
source /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0.2/config/tool_info.txt
source /data5/bsi/bictools/scripts/bior_annotate/v2.4.1//scripts/shared_functions.sh

set -x
VCF=$1
CWD_VCF=`basename $1`
DRILL_FILE=/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0/annotation/drill_file
CATALOG_FILE=/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0/annotation/catalog_file
let count=0

while read drill
do
  let count=$count+1
  SHORT_NAME=`echo "$drill" | cut -f1 `

  CATALOG=`grep -w "$SHORT_NAME" $CATALOG_FILE|cut -f3|head -1`
  CATALOG_COMMAND=`grep -w "$SHORT_NAME" $CATALOG_FILE|cut -f2|head -1`
  TERMS=`echo "$drill" | cut -f2`
  IFS=',' all_terms=( $TERMS )
  separator=" -p "
  drill_opts="$( printf "${separator}%s" "${all_terms[@]}" )"

  if [ -z "$CATALOG" ]
  then
    echo "Error parsing CATALOG. Command used: grep -w \"$SHORT_NAME\" $CATALOG_FILE|cut -f3|head -1"
    exit 100
  fi

  if [ -z "$CATALOG_COMMAND" ]
  then
    echo "Error parsing CATALOG_COMMAND. Command used: grep -w $SHORT_NAME $CATALOG_FILE|cut -f2|head -1"
    exit 100
  fi

  if [ ! -s "$VCF" ]
  then
    echo "$VCF not found. Previous step appears to have failed."
    exit 100
  fi

  cat $VCF | /usr/local/biotools/bior_scripts/2.4.0/bior_pipeline-2.4.0/bin/bior_vcf_to_tjson | /usr/local/biotools/bior_scripts/2.4.0/bior_pipeline-2.4.0/bin/$CATALOG_COMMAND -d $CATALOG | eval /usr/local/biotools/bior_scripts/2.4.0/bior_pipeline-2.4.0/bin/bior_drill ${drill_opts} | /usr/local/biotools/bior_scripts/2.4.0/bior_pipeline-2.4.0/bin/bior_tjson_to_vcf > $CWD_VCF.$count

  START_NUM=`cat $VCF | grep -v '#' | wc -l`
  END_NUM=`cat $CWD_VCF.$count | grep -v '#' | wc -l`
  if [[ ! -s $CWD_VCF.${count} || ! $END_NUM -ge $START_NUM ]]
  then
    /projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/3.0.2//email.sh -f $CWD_VCF.${count} -m annotate.sh -M "bior annotation failed using $CATALOG_FILE" -p $VCF -l $LINENO
    exit 100
  fi

  # Set up for next loop
  if [ "$count" != 1 ]
  then
    let previous=$count-1
    if [[ "NO" == "NO" ]]
    then
      rm $VCF
    fi
  fi

  VCF=${CWD_VCF}.${count}
done <$DRILL_FILE
/usr/bin/perl -pne 's/bior.//g' $CWD_VCF.${count} > ${CWD_VCF}.anno 2>&1 /dev/null

