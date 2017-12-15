### To visualize the states across all chrs for 1 sample 
args<-commandArgs(TRUE)
df<-as.data.frame(read.table(args[1],header=TRUE))
df<-df [(df$chr!="chrX" & df$chr!="chrY"),]
library("GenomicRanges")
# Make GenomicRanges data type from segement.bed file
gr<-makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
gr
source ("http://bioconductor.org/biocLite.R")
#biocLite(c("IRanges","GenomicRanges","ggbio", "biovizBase"))
library("ggbio")
load (system.file("data","hg19IdeogramCyto.rda",package="biovizBase",mustWork=TRUE))
#Get 
data(hg19Ideogram, package="biovizBase")
#seqlengths(gr)<- max(split(end(gr),seqnames(gr)))
seqlengths(gr) <- seqlengths(hg19Ideogram)[names(seqlengths(gr))]
autoplot (seqinfo(gr)[paste0("chr",1:22)])+layout_karyogram(gr,aes(fill=state,color=state))
table(df$state)
