####LOADING NECESSARY LIBRARIES####

library("biomaRt")
library(dplyr)

####PRE-PROCESSING WITH dplyr####

differentialexpression <- read.csv("/home/cilga/diffexpr-results.csv", sep = ",",header = T)
differentialexpression<- differentialexpression[differentialexpression$pvalue > 0.05, ] #pick only transcripts which have p-value more than 0.05
differentialexpression <- differentialexpression[order(-differentialexpression$pvalue),] #order transcriptions by their p-values in descending order.
differentialexpression <- differentialexpression[order(-differentialexpression$log2FoldChange),] #re-order transcripts by their log2foldchange in descending order.

####BIOMART####

ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl") #use ensembl hsapiens dataset for biomaRt.
values<- differentialexpression$Gene #use gene column as values for biomaRt.
data <- getBM(attributes=c("refseq_mrna", "hgnc_symbol"), filters = "refseq_mrna", values = values, mart= ensembl) #do biomaRt between refseq mrna ids and hgnc symbols with using refseq ids as keys.
data$Gene<-data$refseq_mrna 
data$refseq_mrna<-NULL
write.csv(data, file = "analyzed1.csv") #save biomaRt result in a csv file. 


####POST-PROCESSING after BIOMART and merging differential expression data with biomaRt result####

genelistwithhgnc <- merge(differentialexpression,data, by="Gene") #merging biomart results and differential expression data by gene ids.
colnames(genelistwithhgnc)[colnames(genelistwithhgnc)=="data$hgnc_symbol"] <- "hgnc_symbol"  
  
  ##REMOVAL OF UNWANTED COLUMNS WITH nullifying##
  
  genelistwithhgnc$baseMean<-NULL
  genelistwithhgnc$lfcSE<-NULL
  genelistwithhgnc$X<-NULL
  genelistwithhgnc$stat<-NULL
  genelistwithhgnc$LNCaP_10.dsFCS72h_1<-NULL
  genelistwithhgnc$LNCaP_10.dsFCS72h_2<-NULL
  genelistwithhgnc$LNCaP_10.dsFCS72h_3<-NULL
  genelistwithhgnc$X3a_LNCaP_plus_DHT_10.7M_6h_3<-NULL
  genelistwithhgnc$X.home.cilga.Desktop.VPC.rnaseq.tophat_output.SRR1735558.accepted_hits<-NULL
  genelistwithhgnc$X.home.cilga.Desktop.VPC.rnaseq.tophat_output.SRR1735559.accepted_hits<-NULL
  genelistwithhgnc$X.home.cilga.Desktop.VPC.rnaseq.tophat_output.SRR1735560.accepted_hits<-NULL
  genelistwithhgnc$X.home.cilga.Desktop.VPC.rnaseq.tophat_output.SRR1735563.accepted_hits<-NULL
  
genelistwithhgnc <- genelistwithhgnc[order(-genelistwithhgnc$log2FoldChange),] #order merged data frame by descendance of log2foldchange.
mergedtop20<-head(genelistwithhgnc, n = 20, header=T) #pick top 20 rows of merged data frame. 
write.csv(mergedtop20, file = "analyzed1.csv") #save top 20 rows in a csv file. 