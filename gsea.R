y <- read.table("/home/cilga/aa-1.csv",sep = ",",header = T)
head(x)
x$fcsign <- sign(x$log2FoldChange)
x$logP=-log10(x$pvalue)
x$metric= x$logP/x$fcsign
y<-x[,c("Gene", "metric")]
head(y)
y<- cbind(y, x$hgnc_symbol)
y$hgnc_symbol <- y$`x$hgnc_symbol`
y$`x$hgnc_symbol`<-NULL
write.table(y,file="DE_genes12.rnk",quote=F,sep="\t",row.names=F)
x <- read.table("/home/cilga/Desktop/DE_genes12.rnk",sep = "\t",header = T)
y<- read.table("/home/cilga/Desktop/DE_genes11.rnk",sep = "\t",header = T)
