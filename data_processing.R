#### LOADING LIBRARIES ####
library(hash) 
library(dplyr)
library(reshape2)


#### READING METADATA AND MAPPING PATIENT_ID + FILE NAMES ####
metadatareader = read.csv(file = "/home/cilga/Desktop/CILGA/Lack_Lab/TCGA_LEGACY_DATA/TCGA-SARC/metadata-sarc.csv", stringsAsFactors = F, header = T, sep = "\t")
metadatareader$Cancer_Type <- "SARC"

files=hash()    #creating a dictionary with using file name and patient id.
for (i in 1:nrow(metadatareader)){ 
  cancer=metadatareader[i,15]
  type = metadatareader[i,3]
  print(type)
  filename = metadatareader[i,1]
  patientId = paste(c(unlist(strsplit(metadatareader[i,5], split = "-", fixed = T)))[1:4], collapse = "-")
  files[[filename]] = c(metadatareader[i,15], metadatareader[i,3], patientId)
  print(files)
}   
#For loop parses patient id by - and only takes first three seperated part. and writes it in the dictionary. And the nrow means for each row of given file. 
#We need to use unlist because strsplit create a list. and c for combining the vectors and paste for concatenating. 
#Form a list of files with the structural elements like file[filename] = { type, patientId }


#### CREATING A GENE - ISOFORM DICTIONARY ####
genevsisoform <- read.table(file="/home/cilga/Desktop/CILGA/Lack_Lab/TCGA_LEGACY_DATA/TCGA-THYM/genes_result/unc.edu.00f6b730-3424-4a54-9961-551e0ab867ec.2585222.rsem.genes.results", sep = "\t", stringsAsFactors = F, header = T)
#reading gene data from given table in the path. 
genevsisoform$raw_count<- NULL #removal of unwanted columns by equalization them to null.
genevsisoform$scaled_estimate <- NULL

genevsisoformmap=hash() 
#creating dictionary for gene id and isoform ids.
# It loops the first list (into which the geneid is supposed to go)
# then inside, it loops the second list where the isoformId and geneId pair is
# then inside that, it has a check if the id's match, if they do - add the id.
for(i in 1:nrow(genevsisoform)) {  #This loop is required for seperating isoforms from each other.
  isoforms = unlist(strsplit(genevsisoform[i,2], split = ",", fixed = T))
  cat("---------------\n")
  print(isoforms)
  cat("---------------\n")
  for (j in isoforms){ #loop is used for mapping isoforms into a single gene inside of first loop.
    genevsisoformmap[[j]] <- genevsisoform[i,"gene_id"]
    print(j)
  }
}


#### MAPPING ISOFORMS + ISOFORM POOL CREATION ####
isoformfilelistcreator=list.files("/home/cilga/Desktop/CILGA/Lack_Lab/TCGA_LEGACY_DATA/TCGA-SARC/isoforms_result/",full.names = T) #creating a list of files that includes all files in given path. 
isoformfile=read.table("~/Desktop/CILGA/Lack_Lab/TCGA_LEGACY_DATA/TCGA-SARC/isoforms_result/unc.edu.00196819-0436-454b-9a45-5c893f49cd85.2420745.rsem.isoforms.results",sep = "\t",stringsAsFactors = F,header = T)[,1,drop=F]
#read the given file and seperate by tabs and only get first column of the file. 
isoformfile[,"gene_id"] <- "unk" #give unk value for the template file's gene_id column. 

genvsisoformlistKeys = keys(genevsisoformmap) #creating a key for genevsisoform hash map. By this, we will be able to match gene_ids and isofoms list.
unknown_number_tracker=0

#Go over the isoform ids in the template
for (i in 1:nrow(isoformfile)) {        #Go over the isoform ids in gene vs isoform list
  print(isoformfile[i,"isoform_id"])
  if (!(isoformfile[i,"isoform_id"] %in% genvsisoformlistKeys)){   #if isoform cannot be found in keys, it will count the number of unknowns. 
    isoformfile[i,2] <- paste("unknown_gene",unknown_number_tracker,sep = "_")
    unknown_number_tracker = unknown_number_tracker+ 1
  }else{
    isoformfile[i,"gene_id"]=genevsisoformmap[[isoformfile[i,"isoform_id"]]]       
  }
}


 for (i in isoformfilelistcreator){    #create a temporary file for each patient and add their values to the template one in a for loop
  print(i)
  patientspecificisoformfile=read.table(i,sep = "\t",stringsAsFactors = F,header = T)[,2,drop=F] #got the second column of files
  parsedcolnames=paste(c(unlist(strsplit(i, split = ".", fixed = T)))[3], collapse = "-") #create a column named parsedcolnames and parse the filenames by . and get third parsed thing
  parsedcolnames=paste(rev(files[[parsedcolnames]]), collapse = "_") 
  colnames(patientspecificisoformfile)= parsedcolnames
  isoformfile=cbind(isoformfile,patientspecificisoformfile) #make columns of temporary added to the isoformfile
}

#### ADDITION OF TEMPLATE INTO A BIG DATA FRAME ####
TotalTemplate=cbind(cancer, isoformfile)

#### FILTERING DATA FRAME FOR GENE OF INTEREST ####
melted_TotalTemplate <- melt(TotalTemplate) #melting tumor's data frame and drawing boxplot and violin plot.
melted_TotalTemplate <- melted_TotalTemplate[melted_TotalTemplate[,3] == "KDM3B|51780",] #get the given gene id only.ENST00000542866
melted_TotalTemplate$isoformcount<-melted_TotalTemplate$value
melted_TotalTemplate$value<-NULL


#### GENE COUNT POOL CREATION FOR GENE OF INTEREST ####  
totalgenes=list.files("/home/cilga/Desktop/CILGA/Lack_Lab/TCGA_LEGACY_DATA/TCGA-SARC/genes_result/",full.names = T) #creating a list of files that includes all files in given path. 
totalgeneread=read.table("~/Desktop/CILGA/Lack_Lab/TCGA_LEGACY_DATA/TCGA-SARC/genes_result/unc.edu.00196819-0436-454b-9a45-5c893f49cd85.2420744.rsem.genes.results",sep = "\t",stringsAsFactors = F,header = T)[,2,drop=F]

for (i in totalgenes){    #create a temporary file for each patient and add their values to the template one in a for loop
  print(i)
  totalgenereadm=read.table(i,sep = "\t",stringsAsFactors = F,header = T)[,2,drop=F] #got the second column of files
  parsedcolnames=paste(c(unlist(strsplit(i, split = ".", fixed = T)))[3], collapse = "-") #create a column named parsedcolnames and parse the filenames by . and get third parsed thing
  parsedcolnames=paste(rev(files[[parsedcolnames]]), collapse = "*") 
  colnames(totalgenereadm)= parsedcolnames
  totalgeneread=cbind(totalgeneread,totalgenereadm) #make columns of temporary added to the isoformfile
}

totalgeneread$gene_id <-genevsisoform$gene_id
totalgeneread$raw_count<-NULL
totalgeneread$normalized_count<-NULL
melted_genes <- melt(totalgeneread) #melting tumor's data frame and drawing boxplot and violin plot.
melted_genes <- melted_genes[melted_genes[,1] == "KDM3B|51780",] #get the given gene id only.ENST00000542866
melted_genes$gene_id<-NULL
melted_genes$variable<-NULL
melted_genes$totalcount<-melted_genes$value
melted_genes$value<-NULL

repeated_melts <- data.frame()
for (row in 1:nrow(melted_genes)) {
  for (i in 1:3) {
    repeated_melts = rbind(repeated_melts, melted_genes[row, ])
  }
}


#melted_genes <- melt(totalgeneread) #melting tumor's data frame and drawing boxplot and violin plot. #get the given gene id only.
#### ADDITION TO A NEW DATA FRAME ####
genepool=cbind.data.frame(melted_TotalTemplate,repeated_melts)
genepool$totalcount <-genepool$X2004
genepool$X2004 <-NULL


#### MERGING ISOFORM + GENE POOLS ####
FinalDF=data.frame()
FinalDF=rbind(FinalDF,genepool)