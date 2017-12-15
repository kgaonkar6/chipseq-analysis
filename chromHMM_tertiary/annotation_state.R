args<-commandArgs(TRUE)
file<-read.table(args[1],header=T,sep="\t")
overlap_file<-read.table(args[2],header=T,sep="\t")
num_mark_combination<-args[3]


info<-as.data.frame(file[,2:ncol(file)])
rownames(info)<-paste("E",file[,1],sep="")
mean_row<-as.data.frame(rowMeans(info))
rownames(mean_row)<-rownames(info)
min_statenum<-which (mean_row ==min(mean_row))
min_state<-paste("E",min_statenum,sep="")
#final_anno<-c("","")
final_anno_df<-as.data.frame((matrix(ncol=2,nrow=1)))
for (i in 1:nrow(info))
{
state_info<-info[i,]
mark_order<-(colnames(sort(state_info,decreasing=TRUE))[1:num_mark_combination])
#print(paste(mark_order,collapse=":"))
#print(state_info[mark_info,])
state_anno<-c(paste("E",i,sep=""),paste(mark_order,collapse=","))
state_anno_df<-as.data.frame(t(matrix(unlist(state_anno),byrow=T)))
#print (state_anno_df)
final_anno_df<-rbind(final_anno_df,state_anno_df)
}
#final_anno_df

#Get anno for state


overlap<-as.data.frame(overlap_file[1:nrow(overlap_file)-1,2:ncol(overlap_file)])
rownames(overlap)<-paste("E",overlap_file[1:nrow(overlap_file)-1,1],sep="")
new2list<-list('State'=min_state,'Info_Overlap'="Low_Signal")
new2df<-as.data.frame(t(matrix(unlist(new2list),byrow=F)))
overlap<-t(apply(overlap,  2, function(x)(x-min(x))/(max(x)-min(x))))
rownames(overlap)

i<-1
j<-1


for ( j in 1:ncol(overlap) )
{
anno_info<-names(sort(overlap[,j],decreasing=TRUE)[1:3])
values<-paste(anno_info,collapse=",")
state<-colnames(overlap)[j]
anno_values<-c(state,values)



final_anno_df2<-as.data.frame(t(matrix(unlist(anno_values),byrow=T)))
final_anno_df<-rbind(final_anno_df,final_anno_df2)
}


#final_anno_df

#colnames(new2df)<-c("State","Mark_OrderOverlap")
new2df<-unique(final_anno_df)
new2df_agg<-aggregate(V2 ~.,data=new2df,paste,collapse=",")
#as.data.frame(list('State'=new2df_agg$V1,'Info'=new2df_agg$V2))
#colnames(new2df_agg)<-c("State","Mark_Order::Overlap")

for ( i in 1:nrow(new2df_agg))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],",",sep="")
if("RefSeqTSS2kb.hg19.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],",")))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Promoter",sep=",")
}
if(!("enhancer.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],","))) & "K27Ac" %in% unlist(strsplit(new2df_agg$V2[i],",")) & ("K4me3" %in% unlist(strsplit(new2df_agg$V2[i],","))) & "RefSeqTSS2kb.hg19.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],",")))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Active",sep=",")
}
if("enhancer.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],",")) & "K27Ac" %in% unlist(strsplit(new2df_agg$V2[i],",")) & !("Genome.." %in% unlist(strsplit(new2df_agg$V2[i],",")))) 
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Enhancer",sep=",")
}
if( "K27me3" %in% unlist(strsplit(new2df_agg$V2[i],",")) & !("Genome..." %in% unlist(strsplit(new2df_agg$V2[i],","))) & !("RefSeqTSS2kb.hg19.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],","))) & !(("RefSeqTSS.hg19.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],",")))))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Repressed",sep=",")
}
if( "K27me3" %in% unlist(strsplit(new2df_agg$V2[i],",")) & (("K4me3" %in% unlist(strsplit(new2df_agg$V2[i],","))) | ("K4me1" %in% unlist(strsplit(new2df_agg$V2[i],",")))) & "Promoter" %in% unlist(strsplit(new2df_agg$V2[i],",")))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Poissed",sep=",")
}

if("laminB1lads.hg19.bed.gz" %in% unlist(strsplit(new2df_agg$V2[i],",")) & "K9me3" %in% unlist(strsplit(new2df_agg$V2[i],",")) & !("Genome.." %in% unlist(strsplit(new2df_agg$V2[i],","))))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Enhacer",sep=",")
}

if("Genome.." %in% unlist(strsplit(new2df_agg$V2[i],",")))
{
new2df_agg$V2[i]<-paste(new2df_agg$V2[i],"Low Signal",sep=",")
}

new2df_agg$V2[i]<-gsub(",,", ";", new2df_agg$V2[i])
}

final_df<-do.call(rbind.data.frame, Map('c', new2df_agg$V1, new2df_agg$V2))


colnames(final_df)<-c("State","Mark_Order,Overlap_Order;Annotation")
rownames(final_df)<-c(1:nrow(info))
final_df
write.table(final_df,"state_annotation.txt",sep="\t")
