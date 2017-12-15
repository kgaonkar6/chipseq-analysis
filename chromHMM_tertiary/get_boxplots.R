args<-commandArgs(TRUE)
file<-read.table(args[1],header=T)
info<-as.data.frame(file)


library(ggplot2)


png("boxplots.png", units="in", width=2.5, height=2.5, res=100)
par(bg="white")

str(info)
info$state<-factor(info$state,levels=c('E10','E9','E8','E7','E6','E5','E4','E3','E2','E1'),ordered=TRUE)


ggplot(info, aes(y=meanMethyl,x=state))+geom_boxplot()+coord_flip()
dev.off()


