
args<-commandArgs(TRUE)
args[1]

#source("http://bioconductor.org/biocLite.R") 
#biocLite(pheatmap)
library("pheatmap")

#### function

filename <- args[1]
my_data_mat=read.delim(filename, sep="\t", header=F)
### or whatever number you want
#my_data[my_data<=0]<-1
#my_data
#my_data_mat=as.matrix(log2(my_data))
top<-nrow(my_data_mat)
top
genefile<-args[2]
my_gene_dat<-read.delim(genefile,sep="\t",header=F)
rownames(my_data_mat)=paste(my_gene_dat[,1],my_gene_dat[,4], sep="_")
colnames(my_data_mat)=c("1_033","1_037","1_03","1_048","1_052","1_053","1_064","1_119","2_029","2_042","2_045","2_058","2_083","2_087","2_099","2_116","3_043","3_076","AH-IPC","AM-IPC","AO-IPC","B-Tim","D-IPC","foie8")

library("FactoMineR")


res.pca<-PCA (t(my_data_mat), ncp= Inf ,graph=TRUE)
desc_pca<-dimdesc(res.pca)
file_out<-paste("desc_pca",args[3],sep="")
write.infile(desc_pca,file_out,nb.dec=3)
hcpc_xeno<-HCPC(res.pca,nb.clust=-1)
plot(hcpc_xeno,choice="3D.map")
#hcpc_xeno$data.clust
