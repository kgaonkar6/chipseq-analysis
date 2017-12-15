#source("http://bioconductor.org/biocLite.R")
#biocLite("DiffBind")
#biocLite("BioGenerics")

library("DiffBind")
library("locfit")
args<-commandArgs(TRUE)
setwd(args[3])
getwd()
file_name<-paste(args[3],args[4],sep= "/")
file_name
xeno_samples <- read.csv(file_name,header=TRUE)
xeno_samples
xeno <- dba (sampleSheet=xeno_samples)
warnings()
dba.plotPCA(xeno,DBA_CONDITION,label=DBA_CONDITION)
#olap.rate<-dba.overlap(xeno,mode=DBA_OLAP_RATE)
#plot(olap.rate,ylab="# peaks", xlab="Overlap at least this many peaksets")

#need to increase strigency(peak present in 3 or more peakset
mark_output <- paste(args[1],args[2],"corr.pdf", sep= "_" )
path_pdf <- paste (args[3],"report",mark_output, sep= "/")
pdf (path_pdf,width=10,height=10)
xeno <- dba.count (xeno)
dba.plotHeatmap(xeno,correlation=FALSE,scale="row")
#xeno <- dba.contrast (xeno,categories=DBA_CONDITION,minMembers=2)
xeno <- dba.contrast (xeno,group1=1,name1="scram",group2=2,name2="BMI1")
#xeno <- dba.contrast (xeno,group1=1,name1="FRC_K27me3",group2=2,name2="FRD_K27me3")
#xeno <- dba.contrast (xeno,group1=5,name1="FRC2_pol",group2=6,name2="FRD1_pol")

xeno <- dba.analyze(xeno,method=DBA_DESEQ)

#pdf("MAplot.pdf")
#dba.plotMA(xeno)
xeno.db <- dba.report (xeno,DataType=DBA_DATA_FRAME,method=DBA_DESEQ,th=1,contrast=1,bCount=TRUE)
mark_report <- paste("report_",args[1],args[2],".txt", sep= "" )
path_report <-paste(args[3],"report", mark_report, sep="/")
write.table(xeno.db,path_report)
dba.overlap(xeno, xeno$masks$a,mode=DBA_OLAP_RATE)
