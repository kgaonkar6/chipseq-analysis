args<-commandArgs(TRUE)
file<- args[1]
bins<-as.data.frame(read.table(file,header=TRUE))
colnames(bins)
input<-args[2]
input_bins<-as.data.frame(read.table(input,header=TRUE))
colnames(input_bins)


 
# SL1_213_input-3 BIN     SL1_418_input-7 BIN     SL1_485_input-12        BIN     T5_input-7      BIN     T6_input-2      BIN     T7_input-2



#SL1_213_K36me3-7,SL1_418_H3K36me3-6","SL1_485_H3K36me3-11","T5_H3K36me3-6","T6_H3K36me3-1","T7_H3K36me3-1"

library("ggplot2")

#xLabAt=c(0,20,50,80,100)
#xLabel=c("-5Kb","TSS","Genebody","TES","+5Kb")

xLabAt=c(0,20,50,80,100)
xLabel=c("-5Kb","TSS","Genebody","TES","+5Kb")
label=c("T5_H3K36me3","T6_H3K36me3","T7_H3K36me3","SL1_213_K36me3","SL1_418_H3K36me3","SL1_485_H3K36me3")
plot(x=c(1:101),y=(bins$T5_H3K36me3.6-input_bins$T5_input.7),lty=1,type="l",col="blue",ylim=c(0,0.5),axes=FALSE,ylab="IP-INPUT(Reads per million)",xlab="Genomic region",lwd=2)
points(x=c(1:101),y=(bins$T6_H3K36me3.1-input_bins$T6_input.2),lty=1,type="l",col="blue",lwd=2)
points(x=c(1:101),y=(bins$T7_H3K36me3.1-input_bins$T7_input.2),lty=1,type="l",col="blue",lwd=2)
points(x=c(1:101),y=(bins$SL1_418_H3K36me3.6-input_bins$SL1_418_input.7),lty=1,type="l",col="red",lwd=2)
points(x=c(1:101),y=(bins$SL1_213_K36me3.7-input_bins$SL1_213_input.3),lty=1,type="l",col="darkorange",lwd=2)
points(x=c(1:101),y=(bins$SL1_485_H3K36me3.11-input_bins$SL1_485_input.12),lty=1,type="l",col="darkgreen",lwd=2)
#
legend(0,0.5,label,col=c("blue","blue","blue","red","darkorange","darkgreen"), lty=c(1,1,1,1),cex=0.75,lwd=2)



axis(1,at=c(xLabAt),lab=c(xLabel),lwd=1,cex=1, cex.lab=1, cex.axis=1)
axis(2,lwd=1,cex=1, cex.lab=1, cex.axis=1)



