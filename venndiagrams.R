####LOADING LIBRARY####
library(VennDiagram)

####READING PREANALYZED DESEQ2 OUTPUT AND SELECTING SIGNIFICANT ONES####
differentialexpressionswithmart_pre<- read.csv("/home/cilga/Desktop/HT55/DESeq2/preanalyzed/analyzedmmusculus.csv", sep = ",",header = T)
resSig_pre <- differentialexpressionswithmart_pre[ differentialexpressionswithmart_pre$pvalue < .05, ]
resSig_pre<-na.omit(resSig_pre)

down_pre <- resSig_pre[resSig_pre$log2FoldChange < 1, ]
up_pre<- resSig_pre[resSig_pre$log2FoldChange > 1, ]
write.table(down_pre[order(down_pre$log2FoldChange,decreasing=TRUE),], "DownRegulated_pre.txt")
write.table(up_pre[order(up_pre$log2FoldChange,decreasing=TRUE),], "UpRegulated_pre.txt")
down_pre<-na.omit(down_pre)
up_pre<-na.omit(up_pre)

####READING DESEQ2 OUTPUT AND SELECTING SIGNIFICANT ONES####
differentialexpressionswithmart_our<- read.csv("/home/cilga/Desktop/HT55/DESeq2/ourpipeline/analyzedmmusculus.csv", sep = ",",header = T)
resSig_our <- differentialexpressionswithmart_our[ differentialexpressionswithmart_our$pvalue < .05, ]
resSig_our<-na.omit(resSig_our)

up_our <- resSig_our[resSig_our$log2FoldChange > 1, ]
down_our<-resSig_our[resSig_our$log2FoldChange < 1, ]
write.table(up_our[order(up_our$log2FoldChange,decreasing=TRUE),], "UpRegulated_our.txt")
write.table(down_our[order(down_our$log2FoldChange,decreasing=TRUE),], "DownRegulated_our.txt")
down_our<-na.omit(down_our)
up_our<-na.omit(up_our)

####DIAGRAM FOR DOWNREGULATED GENE COUNT####
pdf("downdiagram.pdf")
venn.plot <- venn.diagram(list(down_our$Gene, down_pre$Gene), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Our pipeline", "Preanalyzed"))
grid.draw(venn.plot)
dev.off()

####DIAGRAM FOR UPREGULATED GENE COUNT####
pdf("updiagram.pdf")
venn.plot <- venn.diagram(list(up_our$Gene, up_pre$Gene), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Our pipeline", "Preanalyzed"))
grid.draw(venn.plot)
dev.off()

