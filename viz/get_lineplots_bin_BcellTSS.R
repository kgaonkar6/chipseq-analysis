args<-commandArgs(TRUE)
file<- args[1]
bins<-as.data.frame(read.table(file,header=TRUE))
colnames(bins)
input<-args[2]
input_bins<-as.data.frame(read.table(input,header=TRUE))
colnames(input_bins)


select_bins<-bins[,c("B5_H3K36me3.6","B6_H3K36me3.6","B7_H3K36me3.6","SL1_218_H3K36me3.14","SL1_404_H3K36me3_11","SL1_232_K36me3.6","SL1_269_K36me3_14","SL1_293_K36me3_15","SL1_269_K36me3_14","SL1_73_K36me3.15","SL1_217_K36me3.14")]

#BIN     B5_input-7      BIN     B6_input-7      BIN     B7_input-7      BIN     SL1_218_input-15        BIN     SL1_269_input-15        BIN     SL1_293_input_16        BIN     SL1_404_input_12 BIN     SL1_73_input-1  BIN     SL_217_input-15 BIN     SL_232_input-4


library("ggplot2")

xLabAt=c(0,20,50,80,100)
xLabel=c("-5Kb","TSS","Genebody","TES","+5Kb")
label=c("B5_H3K36me3.6","B6_H3K36me3.6","B7_H3K36me3.6","SL1_218_H3K36me3.14","SL1_404_H3K36me3_11","SL1_232_K36me3.6","SL1_269_K36me3_14","SL1_293_K36me3_15","SL1_269_K36me3_14","SL1_73_K36me3.15","SL1_217_K36me3.14")
plot(x=c(1:101),y=(bins$B5_H3K36me3.6-input_bins$B5_input.7),lty=1,type="l",col="blue",ylim=c(0,0.5),axes=FALSE,ylab="IP-INPUT(Reads per million)",xlab="Genomic region")
points(x=c(1:101),y=(bins$B6_H3K36me3.6-input_bins$B6_input.7),lty=1,type="l",col="blue")
points(x=c(1:101),y=(bins$B7_H3K36me3.6-input_bins$B7_input.7),lty=1,type="l",col="blue")
points(x=c(1:101),y=(bins$SL1_218_H3K36me3.14-input_bins$SL1_218_input.15),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_404_H3K36me3_11-input_bins$SL1_404_input_12),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_232_K36me3.6-input_bins$SL_232_input.4),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_269_K36me3_14-input_bins$SL1_269_input.15),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_293_K36me3_15-input_bins$SL1_293_input_16),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_269_K36me3_14-input_bins$SL1_269_input.15),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_73_K36me3.15-input_bins$SL1_73_input.1),lty=1,type="l",col="red")
points(x=c(1:101),y=(bins$SL1_217_K36me3.14-input_bins$SL_217_input.15),lty=1,type="l",col="red")

legend(0,0.5,label,col=c("blue","blue","blue","red","red","red","red","red","red","red","red"), lty=c(1,1,1,1),cex=0.75,lwd=1)



axis(1,at=c(xLabAt),lab=c(xLabel),lwd=1,cex=1, cex.lab=1, cex.axis=1)
axis(2,lwd=1,cex=1, cex.lab=1, cex.axis=1)



