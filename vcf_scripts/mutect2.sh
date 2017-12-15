tumor=$1
normal=$2
range=$3
sample=$4
/usr/local/biotools/java/jdk1.8.0_20/bin/java -XX:CompileThreshold=1000 -XX:ReservedCodeCacheSize=128m -Xmx10g -Xms5g -jar /data5/bsi/bictools/src/gatk/3.6/GenomeAnalysisTK.jar -T MuTect2 -R $ref -I:tumor $tumor -I:normal $normal -L $range -o $sample
normalSample=`basename $normal | sed -e 's/.bam//g'`
TSample=`basename $tumor | sed -e 's/.bam//g'`
rm $sample.idx
cat $sample | sed -e "s/NORMAL/$normalSample/g" | sed -e "s/TUMOR/$TSample/g" > $sample.tmp

mv $sample.tmp $sample
bgzip -f $sample
tabix -p vcf $sample.gz

