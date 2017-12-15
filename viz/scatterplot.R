args<-commandArgs(TRUE)
args[1]


input<-paste(args[1],"txt",sep=".")
file <- read.table (input,header=T)
filename<-paste(args[1],"FRD7cortest","png",sep=".")
png(filename)
#file<-subset (file, file[,1]<10000)

FRD72<-log10(file[,5]+1.1)
FRC2<-log10(file[,6]+1)
FRC1<-log10(file[,7]+1)
FRD71<-log10(file[,8]+1)

#smoothScatter(FAIRE,DNASE,nbin=20000,main=paste(args[1],"R Square:",round((cor(FAIRE,DNASE,use = "pairwise.complete.obs",method="spearman")^2),3)),cex=3, nrpoints = 0)
smoothScatter(FRD72,FRD71,nbin=1800,cex.axis=2,cex.lab=1.2)
abline(lm(FRD71~FRD72))

filename<-paste(args[1],"FRCcortest","png",sep=".")
png(filename)
smoothScatter(FRC2,FRC1,nbin=1800,cex.axis=2,cex.lab=2)

cor.test(FRD72,FRD71,method="pearson")
cor.test(FRC2,FRC1,method="pearson")
#plot(Dnase,Faire,main=paste(args[1],"R Square:",round((cor(file[,2],file[,1],use = "pairwise.complete.obs")^2),3)),col="blue",pch=20)
abline(lm(FRC1~FRC2)) 

#text(x=c(0.63),y=c(1.09),Corr,3,srt=45,cex=0.75)

