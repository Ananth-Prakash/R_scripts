library(matrixStats)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")

library(mygene)

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Arabidopsis/FROM_PRIDE_FTP")


if ( !exists("EXP_TYPE")) warning("Please specify experiment type variable: EXP_TYE")
## Specify experiment type, i.e. what type of quantification is used, MS1-based or MS2-based. As a rule of thumb, label free and SILAC is MS1-based, and iTRAQ and TMT is MS2-based.
# EXP_TYPE <- "MS2-quant"
EXP_TYPE <- "MS1-quant"


##############
# define FOT normalisation function
# FOT stands for Fraction Of Total. 
# In this normalisation method each protein iBAQ intensity value is scaled to the total amount of signal in a given MS run (column) and transformed to parts per billion (ppb)
fot.normalise <- function(x){
  data.sum <-   apply(x, 2, function(y){sum(y, na.rm=TRUE)})
  # barplot((data.sum), log = "y")
  #
  ##### do ppm normalisation
  x.mat <- as.matrix(x)
  x.mat.ppb <- apply(x.mat, 2, function(i) i/sum(i, na.rm = T) * 1000000000 )
  x.mat.ppb <- as.data.frame(x.mat.ppb)
  colnames(x.mat.ppb) <- paste("ppb.", colnames(x.mat.ppb), sep = "")
  return(x.mat.ppb)
}
##############


##### read  protein expression matrix produced by MaxQuant
tmp  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

## clean up 
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]
## consider those protein groups with more than 1 peptides matched
tmp <- tmp[ tmp[, "Peptides"] > 1, ]

## 
# if experiment is label free use iBAQ quant:
if(EXP_TYPE == "MS1-quant"){
  message("Collecting iBAQ quantification")
  tmp <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]
}


# if experiment is TMT, use reporter intensities, if experiment iTRAQ use Intensities. Note, MaxQuant might change how it reports iTRAQ and TMT intensities.
if(EXP_TYPE == "MS2-quant"){
  message("Collecting MS2 intensities")
  if( any(grepl("Reporter.intensity.corrected", colnames(tmp))) ){
    
    # Filter to consider only protein groups for which the at least 50 % of the sample replicates 
    # have non-zero intensity values
    
    ##### Read sample replicate information
    ##### Note: User defined annotation file
    sample_replicates <- read.table("sample_replicates.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
    
    for(i in 1:nrow(sample_replicates))
    {
      sample_replicates[i,] <- gsub("^","Reporter.intensity.corrected.", sample_replicates[i,], perl=TRUE)
      sample_replicates[i,] <- gsub(" ",".", sample_replicates[i,], perl=TRUE)
      sample_replicates[i,] <- gsub("-",".", sample_replicates[i,], perl=TRUE)
      
      subdata <- tmp[, c(unname(unlist(sample_replicates[i,-1]))), drop=FALSE]
      tmp$replicate_nonzero_instensity_count <- apply(subdata[1:ncol(subdata)], 1, function(x) length(which(x != 0)) )
      
      tmp <- tmp[tmp$replicate_nonzero_instensity_count >= ncol(subdata)/2,]
      
    }
    
    tmp <- tmp[ ,c(2, grep("Reporter.intensity.corrected.[0-9].{1,}.", colnames(tmp))), drop=FALSE ]
  } else {
    tmp <- tmp[ ,c(2, grep("Intensity.", colnames(tmp))), drop=FALSE ]
  }
}


#####
Majority.protein.IDs <- tmp$Majority.protein.IDs
tmp <- tmp[ , -1]

if(EXP_TYPE == "MS1-quant"){
  tmp <- fot.normalise(tmp)  
}
#
tmp <- data.frame( cbind(Majority.protein.IDs, tmp, stringsAsFactors = FALSE) )
##############
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA


data.to.map <- tmp
data.to.map$"ENSG" <- "NA"
data.to.map$"Gene.Symbol" <- "NA"
data.to.map$"unique.gene.count" <- "NA"

##### PXD013868 proteinGroups.txt results taken directly from submitter to PRIDE
##### This has Araport ids.
data.to.map$ENSG <- data.to.map[,c("Majority.protein.IDs")]
data.to.map$ENSG <- gsub("\\.[0-9]","",data.to.map$ENSG, perl=TRUE)
data.to.map$unique.gene.count <- length(unique(strsplit(data.to.map$ENSG, ";")[[1]]))
data.to.map <- data.to.map[data.to.map$unique.gene.count == 1,]
data.to.map$ENSG <- gsub(";.*","",data.to.map$ENSG, perl=TRUE)
data.to.map <- aggregate(data.to.map[ , 2:(ncol(data.to.map)-3)], list("Gene ID" = data.to.map$ENSG), median, na.rm =TRUE)

Ensembl_Biomart_mapping  <- read.table("mart_export.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

merged_data <- merge(x=Ensembl_Biomart_mapping, y=data.to.map,
                     by.x=c("Gene.stable.ID"), by.y=c("Gene ID"),
                     all.x=TRUE, all.y=TRUE)
colnames(merged_data)[1] <- "Gene ID"

write.table(merged_data, "proteinGroups_ppb_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )
##################

