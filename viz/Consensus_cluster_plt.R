args<-commandArgs(TRUE)
args[1]

my_data=read.delim(args[1], sep="\t", header=T)
my_data<-as.matrix(my_data)
library(ConsensusClusterPlus)
filename<-paste(args[2],"Consensus_Cluster.pdf", sep="_")
#pdf(filename)
my_var=apply(my_data,1,var)
my_index=order(my_var, decreasing=T)
nrow(my_data)
top<-0.8*nrow(my_data)
top
my_final_data_mat=my_data[my_index[1:top],]

data_2<-ConsensusClusterPlus(my_final_data_mat, maxK=2, reps=1000, pItem=0.8, pFeature=1,clusterAlg="hc",distance="pearson")


