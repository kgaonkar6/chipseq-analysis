#Script to run multiple correspondence analysis for regions with categorical data

#Rscript mca_mostvariable.R <unionbed from multiple bed files> <top variable regions to study>
#Rplot.pdf will be created visualizing samples in the first 2 most informative components
#Rplot2.pdf will provide a cluster and grouping information
#Rplot4.pdf colors samples according to the grouping in MCA components

#Krutika Gaonkar


args<-commandArgs(TRUE)
#unionbed region with chr<tab>start<tab>stop<tab><per sample(20 sample) data><tab><Count for categorical data>
args[1]
#number of most variable regions
args[2]

#Input file
filename <- args[1]
my_orig_data=read.delim(filename, sep="\t", header=T)

#Reads in categorial data columns
my_data<-my_orig_data[,4:23]
rownames(my_data)=paste(my_orig_data[,1],my_orig_data[,2],my_orig_data[,3], sep="_")

#Reads in count data for categorical data to obtain most variable regions
my_state_data<-my_orig_data[1:3,24:30]
my_data_mat=as.matrix(my_orig_data[,24:30])

#top varaible regions to be studied
top<-args[2]
#calculate variance in count data
my_var=apply(my_data_mat,1,var)
my_index=order(my_var, decreasing=T)
#get top variable regions
my_final_data_mat=my_data[my_index[1:top],]
my_data<-my_final_data_mat


library("FactoMineR")
#create a dataframe with regions as the testing criteria
#samples provide the categorical vales for each region
d<-as.data.frame(t(my_data))
#Run MCA module from Factominer 
res.mca<-MCA(d)
#Obtain statistical summary from mca
summary_mca<-summary.MCA(res.mca,nbelements=Inf,ncp=3,file="summary.out")
#Creates clustering from all the components of the mca 
#Sub groups are created by the module depending on the information retained in each cluster
hcpc_xeno<-HCPC(res.mca,nb.clust=-1)
plot(hcpc_xeno,choice="3D.map")
#Region and cluster per sample is written in file
write.infile(hcpc_xeno$data.clust,"cluster_segments.txt",sep="\t")
