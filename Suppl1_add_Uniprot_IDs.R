# Add uniprot ids to Ensembl gene ids in supplementary file 1.

# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")

library(mygene)

#setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/")
setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/")

#suppl_file <- read.table( "SupplementaryFile_3.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
suppl_file <- read.table( "Gene_distribution_in_organs-GeneNames-Median_bin_values.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


suppl_file$UniProt_ID <- "NA"

#ENSG00000130723, ENSG00000243135 IDs has problems mapping to UniProt ID
filt_file <- suppl_file[suppl_file$GeneID != "ENSG00000130723" & suppl_file$GeneID != "ENSG00000243135",]
filt_entry <- suppl_file[suppl_file$GeneID == "ENSG00000130723" | suppl_file$GeneID == "ENSG00000243135",]

for(i in 1:nrow(filt_file)){
  
  x <- filt_file[i, "GeneID"]
  res <- tryCatch(queryMany(x, scopes="ensembl.gene", fields=c("uniprot"), species=c("human")),
                  error = function(e) {print(0)})
  
  filt_file[ i, "UniProt_ID"] <- paste( unique(unlist(res$uniprot.Swiss.Prot[!is.na(res$uniprot.Swiss.Prot)])), collapse = ";")
  if(filt_file[ i, "UniProt_ID"] == ''){
    filt_file[ i, "UniProt_ID"] <- paste( unique(unlist(res$uniprot.TrEMBL[!is.na(res$uniprot.TrEMBL)])), collapse = ";")
  }
  
  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(filt_file),1)),"%"))
}

filt_file <- rbind(filt_file, filt_entry)
#filt_file <- filt_file[,c(1,2,37,3:36)]
filt_file <- filt_file[,c(1,2,34,3:33)]

filt_file$UniProt_ID <- gsub(";.*","",filt_file$UniProt_ID, perl=TRUE)

write.table(filt_file, file = "SupplementaryFile_3_UniProt_IDs.txt", sep = "\t", row.names = FALSE, quote = FALSE )


filt_file <- read.table( "SupplementaryFile_3_UniProt_IDs.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

ggplot(filt_file, aes(x=present_in_samples))+
  geom_histogram(colour="black", bins=32)+
  theme_bw()+
  xlab("Present in number of organs")+
  ylab("Identified canonical proteins")+
  theme(legend.position="none")+
  ggtitle("Canonical proteins identified in organs")



