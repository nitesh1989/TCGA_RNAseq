#####################################################################
# Nitesh Turaga, nturaga1@jhmi.edu
# Johns Hopkins University

# TCGA- mRNA sequencing data analysis
# Data platform is illumina HiSeq RNAseqV2
# Unpack tar.gz and create expression matrix
#####################################################################

rm(list= ls())
date()
sessionInfo()

#####################################################################
# 1.  Load Libraries
#####################################################################
library(plyr)
library(edgeR)
library(reshape2)


#####################################################################
# 2. Set path
#####################################################################
my.path = "~/Documents/TestRun/TCGA_RNAseqV2" 
setwd(my.path)


#####################################################################
# 3. Unpack directories in this path
#####################################################################

# Set Sample
sample = "COAD"
sample.path = file.path(my.path,"data",sample)

# Add unzipping step
tar.file = list.files(sample.path,pattern = "tar.gz")
untar(tar.file,exdir=sample.path)

# #Read in file_sample_map
# fileMap = read.table("FILE_SAMPLE_MAP.txt",header = T)
# 
# # Add phenotype to file map
# fileMap$sample = substr(x=fileMap$barcode.s.,start =14,stop=15)
# fileMap$phenotype = ifelse(fileMap$sample == "01","Cancer","Normal")
# 
# # files present in level 1 folder
# 
# files = list.files(path = my.path,pattern =".idat",recursive=TRUE,full.names=TRUE)
# files = data.frame(filename = gsub(".+/","",files),path = files)
# 
# # Match the fileMap file with the Level 1 files
# matched.list = plyr::join(x= files ,y=fileMap,by="filename")
# 
# 
# #Make new Cancer and Normal folder
# dir.create("./data")
# dir.create(path="data/Cancer");dir.create(path = "data/Normal")
# 
# #Populate folders
# dataDir = "data"
# 
# # Make matched list character vectors
# matched.list = data.frame(lapply(matched.list,as.character),stringsAsFactors=FALSE)
# 
# 
# # Rename from and to 
# from=matched.list$path
# to = file.path(dataDir,matched.list$phenotype,matched.list$filename)
# 
# 
# # Move files
# file.rename(from= from,to = to)
# 

#####################################################################
# 4. Aggregate Rsem files to a large expression matrix of normalized
# counts
#####################################################################

# Work in sample path
setwd(sample.path)
dir.create("objs")
dir.create("figs")
dir.create("misc")
dir.create("text")

files = dir(full.names=TRUE,recursive=T)

sdrf = files[grep(".sdrf.txt",files)]
files = files[grep(".rsem.genes.normalized_results", files)]


data = lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE,skip=1)

sdrf = read.delim(sdrf,as.is=TRUE)

# Match with sdrf barcodes
barcodes <- sdrf[match(gsub("\\..*$","",gsub(".*unc.edu.","",files)),sdrf[,1]), 2]

# combine all the samples by columns
data <- lapply(1:length(data), function(i) cbind(barcode=barcodes[i],data[[i]]))


# Rename colnames for each list
for(i in 1:length(data)){
    colnames(data[[i]])<-c("barcode","gene_id","normalized_count")
}


# reshape the matrix
rdata <- do.call(rbind, data)
cdata <- dcast(gene_id~barcode,data=rdata, value.var="normalized_count");


filename = paste0("TCGA_RNASeqV2_",sample,"_exprs")
# Save the exprs matrix for the normalized count
save(cdata,file = paste0("objs/",filename,".rda"))

write.csv(cdata, file=paste0("misc/",filename,".csv"),quote=FALSE, row.names=FALSE)

