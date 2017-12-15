###library(gplots)
require(VennDiagram)

K27me3_K4me3 <-as.list(t(read.table ("K27me3_K4me3_genelist.txt", header=F)))
K27me3_K9me3<-as.list(t(read.table ("K27me3_K9me3_genelist.txt", header=F)))
K9me3_K4me3<-as.list(t(read.table ("K9me3_K4me3_genelist.txt", header=F)))
RNAseq<-as.list(t(read.table ("/data2/bsi/tertiary/Urrutia_Raul_rxu01/mrnaseq/xenograft_integration_analysis/list", header=F)))

data_venn<- list(K27_K4me3=K27me3_K4me3,K27_K9me3=K27me3_K9me3,K9_K4me3=K9me3_K4me3,RNAseq=RNAseq)
input<-venn.diagram(list( K27_K4me3=K27me3_K4me3,K27_K9me3=K27me3_K9me3,K9_K4me3=K9me3_K4me3,RNAseq=RNAseq),filename=NULL,fill=2:5,alpha=0.3)
grid.draw(input)

subset (data_venn,K27_K4me3==1 & RNAseq== 1)

