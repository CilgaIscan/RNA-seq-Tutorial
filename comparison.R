####LOADING LIBRARY####
library(compare)

####READING DESEQ2 OUTPUT AND SELECTING UNIQUE ROWS####
ourfulldata <- read.table("/home/cilga/analyzed.csv", sep = ",", header = T)
our <- unique(ourfulldata$hgnc_symbol)

####READING PRE-ANALYZED DESEQ2 OUTPUT AND SELECTING UNIQUE ROWS####
preanalyzedfulldata <- read.table("/home/cilga/preanalyzed.csv", sep = ",", header = T)
pre <- unique(preanalyzedfulldata$hgnc_symbol)

####COMMON GENES BETWEEN TWO ANALYSES####
common <- intersect(pre$hgnc_symbol, our$hgnc_symbol)  
commongenes<-data.frame(common)
unique(commongenes$common)
write.csv(commongenes, "commons.csv")

# subsetour<-subset(ourfulldata, Geneid %in% commongenes$common)
# write.csv(subsetour, "subsettedourdata.csv")
# 
# subsetpre<-subset(preanalyzedfulldata, Gene.ID %in% commongenes$common)
# write.csv(subsetpre, "subsettedpredata.csv")

####UNCOMMON GENES####
uncommon <- pre$hgnc_symbol[!(pre$hgnc_symbol %in% our$hgnc_symbol)]
uncommongenes<-data.frame(uncommon)
unique(uncommon)
write.csv(uncommon, "uncommon.csv")
