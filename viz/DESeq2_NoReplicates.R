
library(DESeq2)
library(Biobase)


dataTable<-read.table("/data2/bsi/staff_analysis/m139467/GeneCount_Westendorf.tsv", sep="\t", header=TRUE, row.names=1)
#head(dataTable)

samples<-data.frame(row.names=c("s_iMAC_CRE_GeneCount","s_iMAC_GFP_GeneCount"), condition=as.factor(c(rep("knockout",1),rep("wild",1))))

#samples<-data.frame(row.names=c("Normal", "Tum"), condition=as.factor(c(rep("Normal",1),rep("Tum",1))))

dds <- DESeqDataSetFromMatrix(countData = dataTable, colData=samples, design=~condition)
#DSeqD_1<-DESeq(dds)

## Simon Anders : running comparisons without replicates ########

rld <- rlogTransformation(dds) #,blind = TRUE)
res <- data.frame(
  assay(rld), 
  avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
  rLogFC = assay(rld)[,2] - assay(rld)[,1] )

my_result<-( res[ order(res$rLogFC), ] )
write.csv(my_result, file="/data2/bsi/staff_analysis/m139467/results_DE_2.txt")
