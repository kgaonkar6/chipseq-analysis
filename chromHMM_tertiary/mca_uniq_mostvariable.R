
args<-commandArgs(TRUE)
#unionbed region
args[1]
#number of most variable regions
args[2]
#source("http://bioconductor.org/biocLite.R") 
#biocLite(pheatmap)
library("pheatmap")

#### function

filename <- args[1]
my_orig_data=read.delim(filename, sep="\t", header=T)
### or whatever number you want
my_data<-my_orig_data[,4:23]
rownames(my_data)=paste(my_orig_data[,1],my_orig_data[,2],my_orig_data[,3], sep="_")


my_state_data<-my_orig_data[1:3,24:30]
my_data_mat=as.matrix(my_orig_data[,24:33])
top<-args[2]

my_var=apply(my_data_mat,1,var)
my_index=order(my_var, decreasing=T)

my_final_data_mat=my_data[my_index[1:top],]
my_data<-my_final_data_mat

#my_data
#my_data=as.matrix(log2(my_data))
#genefile<-args[2]
#my_gene_dat<-read.delim(genefile,sep="\t",header=F)

library("FactoMineR")
d<-as.data.frame(t(my_data))
#is_uniq<-vapply(d,function(x) length(unique(x))>1,logical(1L))
#d<-d[is_uniq]
res.mca<-MCA(d)
#desc_mca<-dimdesc(res.mca,axes=1:3,proba=0.1)

#write.infile(desc_mca,"dimdesc.txt",sep="\t",nb.dec=2)
summary_mca<-summary.MCA(res.mca,nbelements=Inf,ncp=3,file="summary_mca_UNCLiver")
plotellipses(res.mca)
hcpc_xeno<-HCPC(res.mca,nb.clust=-1)
plot(hcpc_xeno,choice="3D.map")
write.infile(hcpc_xeno$data.clust,"cluster_segments.txt",sep="\t")
