
args<-commandArgs(TRUE)
args[1]

#source("http://bioconductor.org/biocLite.R") 
#biocLite(pheatmap)
library("pheatmap")

#### function

filename <- args[1]
my_orig_data=read.delim(filename, sep="\t", header=T)
### or whatever number you want
my_data<-my_orig_data[,4:7]
my_data[my_data<=0]<-1
#my_data
#my_data=as.matrix(log2(my_data))
top<-nrow(my_data)
top
#genefile<-args[2]
#my_gene_dat<-read.delim(genefile,sep="\t",header=F)
rownames(my_data)=paste(my_orig_data[,1],my_orig_data[,2],my_orig_data[,3], sep="_")

library("FactoMineR")


res.pca<-PCA (t(my_data), ncp= Inf ,graph=TRUE)
desc_pca<-dimdesc(res.pca)
plotellipses(res.pca)
write.infile(desc_pca,"desc_pca_K4me1",nb.dec=3)
hcpc_xeno<-HCPC(res.pca,nb.clust=-1)
plot(hcpc_xeno,choice="3D.map")
#hcpc_xeno$data.clust
