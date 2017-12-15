
args<-commandArgs(TRUE)
args[1]

#source("http://bioconductor.org/biocLite.R") 
#biocLite(pheatmap)
library("pheatmap")

library("RColorBrewer")
library("colorRamps")

#### function

filename <- args[1]
plot_color<-blue2red(20)
clu_distance="correlation"
my_data=read.delim(filename, sep="\t", header=F)
### or whatever number you want
my_data[my_data<=0]<-1
#my_data
my_data_mat=as.matrix(log2(my_data[,5:28]))
top<-10000
top
rownames(my_data_mat)=paste(my_data[,1],my_data[,2],my_data[,3],my_data[,4], sep="_")
colnames(my_data_mat)=c("1_033","1_037","1_03_8a","1_048f","1_052","1_053","1_064","1_119","2_029","2_042","2_045","2_058","2_083f","2_087p","2_099","2_116","3_043","3_076","AH-IPC","AM-IPC","AO-IPC","B-Tim","D-IPC","foie8")
colnames(my_data_mat)


my_var=apply(my_data_mat,1,var)
my_index=order(my_var, decreasing=T)

my_final_data_mat=my_data_mat[my_index[1:top],]
colnames(my_final_data_mat)<-colnames(my_data_mat)


my_color=plot_color

basedir <- getwd() 

drows1 <- clu_distance
dcols1 <- clu_distance

mark=args[2]

filename <- paste(mark, "_top",top, "_heatmap.2way.pdf", sep="")
outfile <- paste(basedir, filename, sep="/")



hm.parameters <- list(my_final_data_mat, 
  color = plot_color,
  scale = "row",
  kmeans_k = NA,
  cutree_rows=3,	
  show_rownames = F, show_colnames = T,
  main = " ",
  clustering_method = "complete",
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = drows1, 
  clustering_distance_cols = dcols1)


do.call("pheatmap", c(hm.parameters, filename=outfile))

mark_tree<- do.call("pheatmap", c(hm.parameters))
mark.clust<-cbind(my_final_data_mat,cluster=cutree(mark_tree$tree_row, k=3))
mark.clust


