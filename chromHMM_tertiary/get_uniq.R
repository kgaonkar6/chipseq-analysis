args<-commandArgs(TRUE)
args[1]

#source("http://bioconductor.org/biocLite.R")
#biocLite(pheatmap)
library("pheatmap")

#### function

filename <- args[1]
my_orig_data=read.delim(filename, sep="\t", header=T)
### or whatever number you want
my_data<-my_orig_data[,4:24]
#my_data[my_data!=2]<-"N"
#my_data[my_data==2]<-"Y"

#my_data
#my_data=as.matrix(log2(my_data))
top<-nrow(my_data)
top
#genefile<-args[2]
#my_gene_dat<-read.delim(genefile,sep="\t",header=F)
index_row<-c(1:top)
rownames(my_data)=paste(my_orig_data[,1],my_orig_data[,2],my_orig_data[,3], sep="_")

library("FactoMineR")
d<-as.data.frame(t(my_data))
d
#is_uniq<-vapply(d,function(x) length(unique(x))>2,logical(1L))
#d<-d[is_uniq]
#dim(d)
