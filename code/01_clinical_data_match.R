#####################################################################
# Nitesh Turaga, nturaga1@jhmi.edu
# Johns Hopkins University

# TCGA- mRNA sequencing data analysis
# Data platform : illumina HiSeq RNAseqV2
# Match with Clinical Data

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
# 3. Sample working directory
#####################################################################

# Set Sample
sample = "COAD"
sample.path = file.path(my.path,"data",sample)

load("objs/TCGA_RNASeqV2_COAD_exprs.rda")

dim(cdata)


# Make rownames the gene_id
rownames(cdata) = cdata$gene_id
cdata$gene_id = NULL

# Make the column names contain Tumor vs Normal information
colnames(cdata) = substr(colnames(cdata),1,16)

# Add phenotype to file map
phenotype = data.frame(Sample = colnames(cdata),status = substr(x=colnames(cdata),start =14,stop=15))

# Save expression set as object

save(cdata,file = paste0("objs/TCGA_RNASeqV2_",sample,"_exprs2.rda"))
save(phenotype,file = paste0("objs/",sample,"_tumor_normal_status.rda"))

# Load the clinical patient data
patient.data = list.files(sample.path,pattern = "clinical_patient")
patient.data = read.table(patient.data,fill=T,stringsAsFactors=F,sep ="\t",header =T)
patient.data = patient.data[c(-1,-2),]

# Match the data
subset.tumor = patient.data[which(patient.data$bcr_patient_barcode %in% substr(colnames(cdata),1,12)),]


#Load normal patient data
normal.data = list.files(sample.path,pattern = "normal_control")
normal.data = read.table(normal.data,fill=T,stringsAsFactors=F,sep ="\t",header =T)
normal.data = normal.data[c(-1),]

subset.normal = normal.data[which(normal.data$bcr_patient_barcode %in% substr(colnames(cdata),1,12)),]



dim(subset.data)

table(subset.data$gender)
table(subset.data$race)
table(subset.data$ethnicity)
table(subset.data$tissue_source_site)

write.csv(subset.data,file = "objs/_expression_matched.csv")



    


