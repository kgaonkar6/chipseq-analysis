
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
rownames(my_data)=paste(my_orig_data[,1],my_orig_data[,2],my_orig_data[,3], sep="_")

library("FactoMineR")
d<-as.data.frame(t(my_data))
is_uniq<-vapply(d,function(x) length(unique(x))>1,logical(1L))
d<-d[is_uniq]
d
res.mca<-MCA(d,ncp=3)
desc_mca<-dimdesc(res.mca,axes=1:2,prob=1)
#library("factoextra")
#fviz_screeplot(res.mca,ncp=10)
write.infile(desc_mca,"dimdesc.txt",sep="\t",nb.dec=2)
#summary_mca<-summary.MCA(res.mca,nbelements=Inf,ncp=3,file="summary_mca_UNCLiver")
plotellipses(res.mca)
hcpc_xeno<-HCPC(res.mca,nb.clust=-1)
plot(hcpc_xeno,choice="3D.map")
write.infile(hcpc_xeno$data.clust,"cluster_segments.txt",sep="\t",nd.dec=2)
