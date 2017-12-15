## get_pvclust.R [pvclust input file] [mark] [pvclust output file] ##

args<-commandArgs(TRUE)

#The scripts take a numberic matrix file to create pvclust plot between the control and IP

library(pvclust)


#pvclust_input=argv[1]
args[1]

#mark=argv[2]
args[2]
mark <- args[2]

#pvclust_output=argv[3]
args[3]
filename <-args[3]



#name<-unlist(strsplit (args[1],"[.]"))
#filename <-paste(name[1],"pdf",sep=".")
#mark<-unlist(strsplit (name, "[_]"))

pdf (filename,width=10,height=7.5)
file <- read.table(args[1], header=TRUE)
par(bg = 'white')
#results <- pvclust ( file , method.dist="cor", method.hclust="average", nboot=1000)
results <- pvclust (file)
plot(results,cex=0.75,main=mark)


