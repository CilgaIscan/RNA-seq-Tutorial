#### LOADING LIBRARIES ####
library(ggplot2) 
library(dplyr)
library(reshape2)
library(grid)

FinalDF<- allinoneend

#### SELECTING NECESSARY COLUMNS + ISOFORM/TOTALCOUNT VALUE CALCULATION ####
selection=select(FinalDF, contains("cancer"), contains("gene_id"), contains("isoform_id"), contains ("variable"), contains("_Primary Tumor"), contains("isoformcount"), contains("totalcount"))
selection$value=selection$isoformcount/selection$totalcount


#### PLOTTING ####
plot.new()

melted_graph <- rbind.data.frame(FinalDF) #binding rows of melted selection and normal in a single frame called as melted_total
melted_graph$sample <- melted_graph$variable #and we renamed variable column as sample
melted_graph$sample <- as.character(melted_graph$sample) #we claimed that sample column has characters
melted_graph <- melted_graph[melted_graph[2] == "uc003lcy.1",]
melted_graph$sample=unlist(lapply(X = melted_graph[,"sample"], 
                                  FUN = function(x) unlist(strsplit(x, fixed=T,split='_'))[2])) 
melted_graph <- melted_graph[melted_graph[7] == "Primary Tumor",]

A=ggplot(data = melted_graph, aes(x=cancer,y=isoformcount/totalcount, colour=cancer,fill=cancer)) + geom_violin() + geom_boxplot(width=0.1,notch = T, color="Black")
#+facet_wrap(~ isoform_id)
A1=A+ggtitle("TCGA Tumor Isoform One")
A2=A1+coord_flip()
ggsave("TCGA_legacy_tumor_i1.pdf", scale = 1, width = 11, height = 8, units = "in")
