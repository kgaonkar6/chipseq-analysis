library(ChIPpeakAnno)
hg19.genes.bed<-read.delim("/data2/external_data/Kocher_Jean-Pierre_m026645/s200625.lymphoma_SL1/anno/hg19.iGenome.Genes.bed", header=F)
names (hg19.genes.bed)<-c("chr","start","end","gene_id","strand")
annoData.refGenes<-makeGRangesFromDataFrame(hg19.genes.bed,keep.extra.columns=T)
names (annoData.refGenes)<-hg19.genes.bed$gene_id

union_seg_df<-as.data.frame(read.table("union_segmentation_10_uniq_no_sexChr.bed", header=T, sep="\t"))

union_seg<-makeGRangesFromDataFrame(union_seg_df,keep.extra.columns=T)


anno.refgenes<-annotatePeakInBatch(union_seg, AnnotationData=annoData.refGenes, output="both", ignore.strand=F)

write.table(as.data.frame(anno.refgenes),file="Bcell_control_patient_10.txt")


