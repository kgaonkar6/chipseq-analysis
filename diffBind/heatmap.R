##Usage ./heatmap [down/up] <numeric table for all the groups> <pdf name> ##

args<-commandArgs(TRUE)
args[1]
args[2]
args[3]

file <- read.table (args[1],header=F,sep="\t")

dat <<- as.matrix(file)
library (gplots)
#rownames(dat) <- read.table (args[2],sep="\t")

#rownames(dat) <- c("H3K4me3","H3K4me1","H3K27ac","Foxo3a","H3K9me3","H3K27me3")
dat
num<-0
name<-unlist(strsplit (args[1],"[.]"))


filename <- args[3]
pdf(filename)
par(bg = 'white')
hm<-heatmap.2 (dat, Colv=NA,Rowv=NA, col=greenred(300),trace="none",labRow="")
as.hclust(hm$rowDendrogram)

