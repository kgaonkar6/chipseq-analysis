args<-commandArgs(TRUE)
file<-read.table(args[1],header=TRUE)
info<-as.data.frame(file)


library(ggplot2)


png("boxplots.png", units="in", width=11, height=8.5, res=300)
par(bg="white")
qplot( name,log10(as.numeric(info$Cov)),data=info, geom="boxplot" )
dev.off()


