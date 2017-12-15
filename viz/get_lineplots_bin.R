args<-commandArgs(TRUE)
file<- args[1]
bins<-as.data.frame(read.table(file,header=TRUE))
colnames(bins)
input<-args[2]
input_bins<-as.data.frame(read.table(input,header=TRUE))
colnames(input_bins)


 
#BIN     B5_input-7      BIN     B6_input-7      BIN     B7_input-7      BIN     SL1_218_input-15        BIN     SL1_269_input-15        BIN     SL1_293_input_16        BIN     SL1_404_input_12 BIN     SL1_73_input-1  BIN     SL_217_input-15 BIN     SL_232_input-4


#B5__H3K27me3-5","B6__H3K27me3-5","B7__H3K27me3-5","SL1_218_H3K27me3-13","SL1_232_K27me3-5","SL1_269_K27me3_13","SL1_293_K27me3_14","SL1_404_H3K27me3_10","SL1_73_K27me3-14","SL1_217_K27me3-13"

library("ggplot2")

#xLabAt=c(0,20,50,80,100)
#xLabel=c("-5Kb","TSS","Genebody","TES","+5Kb")

xLabAt=c(0,50,100)
xLabel=c("-5KB","TSS","+5KB")
label=c("B5_H3K27me3","B6_H3K27me3","B7_H3K27me3","SL1_218_H3K27me3","SL1_404_H3K27me3","SL1_232_K27me3","SL1_269_K27me3","SL1_293_K27me3","SL1_73_K27me3","SL1_217_K27me3")
plot(x=c(1:101),y=(bins$B5__H3K27me3.5-input_bins$B5_input.7),lty=1,type="l",col="blue",ylim=c(0,0.5),axes=FALSE,ylab="IP-INPUT(Reads per million)",xlab="Genomic region",lwd=2)
points(x=c(1:101),y=(bins$B6__H3K27me3.5-input_bins$B6_input.7),lty=1,type="l",col="blue",lwd=2)
points(x=c(1:101),y=(bins$B7__H3K27me3.5-input_bins$B7_input.7),lty=1,type="l",col="blue",lwd=2)
points(x=c(1:101),y=(bins$SL1_218_H3K27me3.13-input_bins$SL1_218_input.15),lty=1,type="l",col="red",lwd=2)
points(x=c(1:101),y=(bins$SL1_404_H3K27me3_10-input_bins$SL1_404_input_12),lty=1,type="l",col="darkorange",lwd=2)
points(x=c(1:101),y=(bins$SL1_232_K27me3.5-input_bins$SL_232_input.4),lty=1,type="l",col="darkgreen",lwd=2)
points(x=c(1:101),y=(bins$SL1_269_K27me3_13-input_bins$SL1_269_input.15),lty=1,type="l",col="darkviolet",lwd=2)
points(x=c(1:101),y=(bins$SL1_293_K27me3_14-input_bins$SL1_293_input_16),lty=1,type="l",col="gold4",lwd=2)
points(x=c(1:101),y=(bins$SL1_73_K27me3.14-input_bins$SL1_73_input.1),lty=1,type="l",col="coral",lwd=2)
points(x=c(1:101),y=(bins$SL1_217_K27me3.13-input_bins$SL_217_input.15),lty=1,type="l",col="deeppink",lwd=2)

legend(0,0.5,label,col=c("blue","blue","blue","red","darkorange","darkgreen","darkviolet","gold4","coral","deeppink"), lty=c(1,1,1,1),cex=0.75,lwd=2)



axis(1,at=c(xLabAt),lab=c(xLabel),lwd=1,cex=1, cex.lab=1, cex.axis=1)
axis(2,lwd=1,cex=1, cex.lab=1, cex.axis=1)



