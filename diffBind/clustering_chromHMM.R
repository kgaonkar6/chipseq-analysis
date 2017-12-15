
args<-commandArgs(TRUE)
args[1]

#source("http://bioconductor.org/biocLite.R") 
#biocLite(pheatmap)
library("pheatmap")

library("RColorBrewer")
library("colorRamps")


#### function

#clu_distance="correlation"
plot_color=blue2red(7)

my_data=read.delim(args[1], sep="\t", header=T)
my_data
my_data_mat=as.matrix(my_data[,4:27])
my_data_mat
rownames(my_data_mat)=paste(my_data[,1],my_data[,2],my_data[,3], sep="_")

my_var=apply(my_data_mat,1,var)
my_var
my_index=order(my_var, decreasing=T)
ind <- apply(my_data_mat, 1, var) == 0
my_data_mat <- my_data_mat[!ind,]

#my_var=apply(my_data_mat,1,var)
#my_var
#my_index=order(my_var, decreasing=T)
#my_final_data_mat=as.numeric(my_data_mat[my_index[1:10000],])


#my_final_data_mat=as.numeric(my_data_mat)
#my_final_data_mat

#heatmap.2(my_data_mat)
my_color=plot_color
basedir <- getwd() 

clu_distance<-dist(my_data_mat)
drows1 <- clu_distance
dcols1 <- clu_distance

mark=args[2]

filename <- paste(mark, "_heatmap.2way.pdf", sep="")
outfile <- paste(basedir, filename, sep="/")

bk<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
plot_color<-c("red","pink","magenta","coral","brown","green","forestgreen","aquamarine","orange","yellow","khaki4","cyan","deepskyblue","gray","blue")

hm.parameters <- list(my_data_mat, 
  color = plot_color,
#  scale = "row",
  kmeans_k = NA,breaks=NA,
  show_rownames = F, show_colnames = T,
  main = " ",
  clustering_method = "complete",
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = drows1, 
  clustering_distance_cols = dcols1,
   na.rm=T)


do.call("pheatmap", c(hm.parameters, filename=outfile))

#}



#### process data

#my_file=args[1]

#get_heatmap(filename=my_file, top=20000, clu_distance="correlation", plot_color=blue2red(100))

