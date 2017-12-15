
args<-commandArgs(TRUE)
args[1]

#source("http://bioconductor.org/biocLite.R") 
#biocLite(pheatmap)
library("pheatmap")

library("RColorBrewer")
library("colorRamps")

#### function

filename <- args[1]
plot_color<-c("red","orange","pink","brown","blue","green","gray")
clu_distance="correlation"
my_data=read.delim(filename, sep="\t", header=T)
### or whatever number you want
my_data[my_data<=0]<-1
#my_data
my_data_mat=as.matrix(my_data[,9:28])
top<-nrow(my_data_mat)
top

rownames(my_data_mat)=paste(my_data[,1],my_data[,2],my_data[,3], sep="_")
colnames(my_data_mat)=c("366","389","392","393","394","396","398","399","400","402","406","408","409","410","412","413","414","416","417","420")

my_var=apply(my_data_mat,1,var)
my_index=order(my_var, decreasing=T)


my_final_data_mat=my_data_mat[my_index[1:top],]
colnames(my_final_data_mat)<-colnames(my_data_mat)
cancer_type<-data.frame(factor(c("Cirrhosis","Control","Cirrhosis","Control","Cirrhosis","Control","Cirrhosis","Cirrhosis","Control","Cirrhosis","SevereAH","SevereAH","SevereAH","SevereAH","SevereAH","SevereAH","SevereAH","SevereAH","Cirrhosis","Cirrhosis")))
colnames(cancer_type)<-"Type"
rownames(cancer_type)<-colnames(my_final_data_mat)
cancer_type
cancer_type_color<-list(ID=c(SevereAH="red90",Cirrhosis="navy",Control="gray"))

my_color=plot_color







basedir <- getwd() 

drows1 <- clu_distance
dcols1 <- clu_distance

mark=args[2]

filename <- paste(mark, "_top",top, "_heatmap.2way.pdf", sep="")
outfile <- paste(basedir, filename, sep="/")

quantile.range<- quantile(my_final_data_mat,probs=seq(0,1,0.01))
bk<- seq(1:8)

plot_color<-blue2red(length(bk)-1)

hm.parameters <- list(my_final_data_mat,  
  color=my_color, 
  breaks=bk,
  annotation_col=cancer_type,
  annotation_colors=cancer_type_color[1],
  annotation_legend=TRUE,
  kmeans_k = NA,
  cutree_rows=3,	
  show_rownames = F, show_colnames = T,
  scale="none",
  dendrogram="row",
  main = mark,
  cluster_rows = T, cluster_cols = T)


do.call("pheatmap", c(hm.parameters, filename=outfile))


mark_tree<- do.call("pheatmap", c(hm.parameters))
mark.clust<-data.frame(cutree(mark_tree$tree_row, k=3))
rownames (mark.clust)<-rownames(my_final_data_mat)
filename <- paste(mark, "_top",top, "_cluster.txt", sep="")
outfile <- paste(basedir, filename, sep="/")
write.table(mark.clust,file=outfile)


