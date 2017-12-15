args<-commandArgs(TRUE)
args[1]

library("pheatmap")
library("RColorBrewer")
library("colorRamps")



#### function

get_heatmap=function(filename=my_file, top=20000, clu_distance="correlation", plot_color=blue2red(100)) {

my_data=read.delim(filename, sep="\t", header=T)

my_data_mat=as.matrix(log2(my_data[,5:28]+1))
rownames(my_data_mat)=paste(my_data[,1],my_data[,2],my_data[,3],my_data[,4], sep="_")

my_var=apply(my_data_mat,1,var)
my_index=order(my_var, decreasing=T)

my_final_data_mat=my_data_mat[my_index[1:top],]

my_color=plot_color

basedir <- getwd() 

drows1 <- clu_distance
dcols1 <- clu_distance

sample_name=args[2]

filename <- paste(sample_name, "_top",top, "_heatmap.2way.pdf", sep="")
outfile <- paste(basedir, filename, sep="/")


hm.parameters <- list(my_final_data_mat, 
  color = plot_color,
  scale = "row",
  kmeans_k = NA,
  show_rownames = F, show_colnames = T,
  main = " ",
  clustering_method = "complete",
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = drows1, 
  clustering_distance_cols = dcols1)


do.call("pheatmap", c(hm.parameters, filename=outfile))

}



#### process data

my_file=args[1]

get_heatmap(filename=my_file, top=20000, clu_distance="correlation", plot_color=blue2red(100))

