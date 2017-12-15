##Usage ./heatmap [down/up] <numeric table for all the groups> <pdf name> ##

args<-commandArgs(TRUE)
args[1]
args[2]
args[3]

file <- read.table (args[1],header=F,sep="\t")

dat <<- as.matrix(file)
dat

library (gplots)
rownames(dat) <- as.matrix(read.table (args[2]),sep="\t",header=F)
rownames(dat)

colnames(dat) <- c("CP5","CP25","KO25")
dat
num<-0


name <- args[3]
filename<-paste (name,"pdf",sep=".")
pdf(filename,width=15, height=15)
par(bg = 'white')
hm<-heatmap.2 (dat, Colv=NA,Rowv=NA, col=greenred(300),trace="none",scale="row",cexCol=1.3,main=name, cexRow=0.7,lhei=c(0.4,5),lwid=c(1.5,2.5))
as.hclust(hm$rowDendrogram)

