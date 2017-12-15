# Simple Pie Chart
slices <- c(3,12,11,22,19,3,4,112)
lbls <- c("ProximalPromoter","Promoter1k","Promoter3k", "Upstream10k","Downstream10k","Downstream1k","Downstream3k","Genebody")
pie(slices, labels = lbls, main="Genome wide distribution of binding")


#"Proximal" 3
#"Promoter1k" 12
#"Promoter3k" 11
#"Upstream10k" 22
#"Downstream10k" 19
#"Downstream1k" 3
#"Downstream3k" 4
#"Genebody" 112





#window 1000

#Promoter1k	39
#ProximalPromoter	14
#Promoter3k	94
#Genebody	694
#Pericentromere	9
#OtherIntergenic	585
#Genedesert	109


#Promoter1k	4048
#ProximalPromoter	1436
#Promoter3k	6705
#Genebody	51693
#Pericentromere	11
#OtherIntergenic	29999
#Genedesert	3153

#Window 200

#Promoter1k	5
#ProximalPromoter	2
#Promoter3k	44
#Genebody	87
#Pericentromere	12
#OtherIntergenic	175
#Genedesert	19


#Promoter1k	54
#Promoter3k	60
#ProximalPromoter	10
#Genebody	270
#Pericentromere	3
#OtherIntergenic	158
#Genedesert	7


