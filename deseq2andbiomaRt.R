#### LOADING LIBRARIES
library(DESeq2)
library(RColorBrewer)
library(gplots)

#### LOADING DATASET ####
CountTable <- read.table("/home/cilga/Desktop/mmusculus/mmusculus_featurecounts/rnaseq_mm.txt",sep = "\t",header=TRUE, fill = T)
rownames(CountTable) <- CountTable[,1]
CountTable <- CountTable[,-1]
CountTable <- CountTable[ ,6:ncol(CountTable)]
CountTable$FDR<-NULL

## Converting data frame into matrix ##
CountTable <- as.matrix(CountTable)
head(CountTable)

## Assign conditions## [rep("name", times)] ##
conditions <- factor(c(rep("female", 2), rep("male", 2)), levels=c("female","male"))


#### DESEQ2 ####
## Create a coldata frame and instantiate the DESeqDataSet ##
(coldata <- data.frame(row.names=colnames(CountTable), conditions))
datasetforDESEQ <- DESeqDataSetFromMatrix(countData=round(CountTable), colData=coldata, design=~conditions)
datasetforDESEQ

## Run the DESeq pipeline ##
datasetforDESEQ <- DESeq(datasetforDESEQ)


####PLOTTING####
## Plot dispersions ##
png("qc-dispersions-mm-ourpipeline.png", 1000, 1000, pointsize=20)
plotDispEsts(datasetforDESEQ, main="Dispersion plot")
dev.off()

## Regularized log transformation for clustering/heatmaps, etc ##
rld <- rlogTransformation(datasetforDESEQ)
head(assay(rld))
hist(assay(rld))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(conditions))])

#Colors for plots below
#Ugly:
#(mycols <- 1:length(unique(condition)))

## Sample distance heatmap ##
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples-mm-ourpipeline.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[conditions], RowSideColors=mycols[conditions],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

## Principal components analysis ##
# Could do with built-in DESeq2 function : DESeq2::plotPCA(rld, intgroup="condition")
rld_pca <- function (rld, intgroup = "conditions", ntop = 500, colors=NULL,
                     legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste,
                     collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }
    else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab,
       ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)),
                                    cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
           #pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
          #terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca-mm-ourpipeline.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="conditions", xlim=c(-75, 35))
dev.off()

# Get differential expression results ##
res <- results(datasetforDESEQ)
table(res$padj<0.05)

## Order by adjusted p-value ##
res <- res[order(res$padj), ]

# Merge with normalized count data ##
resdata <- merge(as.data.frame(res), as.data.frame(counts(datasetforDESEQ, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


#### WRITE RESULTS ####
## Write results ##
write.csv(resdata, file="diffexpr-results-mm-ourpipeline.csv")

## Examine plot of p-values ##
hist(res$pvalue, breaks=50, col="grey")


####BIOMART#####
####LOADING NECESSARY LIBRARIES####
library("biomaRt")
library(dplyr)
library(calibrate)


####PRE-PROCESSING WITH dplyr####
differentialexpression <- read.csv("/home/cilga/diffexpr-results-mm-ourpipeline.csv", sep = ",",header = T)


####BIOMART####
ensembl=useMart("ensembl")
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
values<- differentialexpression$Gene #use gene column as values for biomaRt.
data <- getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", values = values, mart= ensembl) #do biomaRt between refseq mrna ids and hgnc symbolswith using refseq ids as keys.
data$Gene<-data$ensembl_transcript_id
data$ensembl_transcript_id<-NULL
write.csv(data, file = "analyzed.csv") #save biomaRt result in a csv file.


####POST-PROCESSING after BIOMART and merging differential expression data with biomaRt result####
genelistwithhgnc <- merge(differentialexpression,data, by="Gene") #merging biomartresults and differential expression data by gene ids.
colnames(genelistwithhgnc)[colnames(genelistwithhgnc)=="data$gene_id"] <"hgnc_symbol"
genelistwithhgnc <- genelistwithhgnc[order(genelistwithhgnc$log2FoldChange),] #ordermerged data frame by descendance of log2foldchange.
write.csv(genelistwithhgnc, file = "analyzedmmusculus.csv") #save top 20 rows in acsv file.
mergedtop20<-head(genelistwithhgnc, n = 20, header=T) #pick top 20 rows of mergeddata frame.
write.csv(mergedtop20, file = "analyzedmmusculus-top20.csv") #save top 20 rows in a csv file.


#### READING DESEQ2 OUTPUT WITH HGNC GENE SYMBOLS ####
res <- read.table(file = "/home/cilga/analyzedmmusculus.csv", header=TRUE, sep = ",")
head(res)
png("diffexpr-volcanoplot-mm-ourpipeline.png", 1200, 1000, pointsize=20)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20,col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue),pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=hgnc_symbol, cex=.8))
dev.off()

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=hgnc_symbol,
                                          cex=textcx, col=2))
  }
}
png("diffexpr-maplot-mm-ourpipeline.png", 1500, 1000, pointsize=20)
maplot(res, main="MA Plot")
dev.off()
