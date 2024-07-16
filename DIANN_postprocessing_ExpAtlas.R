#### DIA-NN postprocessing
### For files to Expression Atlas

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd('/Users/ananth/Documents/DIANN/')

DIA_datasets <- data.frame(Dataset = c(
  "PXD012254",
  "PXD001506",
  "PXD001764",
  "PXD025705",
  "PXD004684",
  "PXD032076",
  "PXD025431",
  "PXD002732",
  "PXD000672",
  "PXD004873",
  "PXD019594",
  "PXD022872",
  "PXD034908",
  "PXD018830",
  "PXD039665",
  "PXD033060",
  "PXD031419"
))

for(i in 1:nrow(DIA_datasets )){
  Dataset <- DIA_datasets[i,]
  #Dataset <- "PXD031419"

  tmp  <- read.table(file=paste("/Users/ananth/Documents/DIANN/",Dataset,"/SwissProt/report_AllSamples.tsv",sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  tmp$File.Name <- gsub(".*PXD","PXD",tmp$File.Name,perl=TRUE)
  tmp$File.Name <- gsub("/submitted.*","",tmp$File.Name,perl=TRUE)
  tmp <- tmp[tmp$File.Name != "",]
  tmp <- tmp[!grepl("Library|POOL",tmp$Run, ignore.case=TRUE),]
  # For dataset PXD031419
  tmp <- tmp[!grepl("20161118_RiMod_H_10.|20161121_RiMod_I_10.|20161121_RiMod_I_5.",tmp$Run, ignore.case=TRUE),]
  

  #tmp_submatrix <- tmp[,c("File.Name","Run","Genes","Genes.Normalised","Modified.Sequence","Stripped.Sequence","Precursor.Id")]
  tmp_submatrix <- tmp[,c("File.Name","Run","Genes","Genes.Normalised","Stripped.Sequence")]
  tmp_submatrix$Genes.Normalised <- as.numeric(tmp_submatrix$Genes.Normalised)
  tmp_submatrix <- unique(tmp_submatrix)
  tmp_submatrix <- spread(tmp_submatrix, Run, Genes.Normalised)

  # don't consider peptides which map to more than one gene (gene groups)
  # peptide mappings should be unique to a gene
  tmp_submatrix <- tmp_submatrix[!grepl(";",tmp_submatrix$Genes),]
  #Remove contaminants
  tmp_submatrix <- tmp_submatrix[!grepl("SWISS-PROT",tmp_submatrix$Genes),]
  tmp_submatrix <- tmp_submatrix[tmp_submatrix$Genes!= "" & tmp_submatrix$Genes!="(Bos" & tmp_submatrix$Genes!="(S.avidinii)" & tmp_submatrix$File.Name!="",]
  
  PSMs_per_gene <- tmp_submatrix  %>%
   group_by(File.Name,Genes) %>%  #Protein.Group can be removed to give accurate count
   summarise(Total_PSMs = n(),
            Distinct_StrippedSeq_PSMs = n_distinct(Stripped.Sequence))

  # Check if to filter by "Total peptides" or "Distinct peptides"?
  # Total_PSMs  = Total number of peptidoforms mapped to a gene (can have duplicates)
  # Distinct_PSMs = Total number of distinct peptidoforms only that are mapped to a gene (does not have duplicates) 

  ### Select only those genes which have evidences of more than 1 PSM
  #Filtered_PSMs_per_gene <- PSMs_per_gene[PSMs_per_gene$Total_PSMs > 1,]
  Filtered_PSMs_per_gene <- PSMs_per_gene[PSMs_per_gene$Distinct_StrippedSeq_PSMs > 1,]

  #all_genes <- data.frame(Filtered_PSMs_per_gene$Genes)
  #Filtered_genes <- data.frame(Genes=unique(all_genes$Filtered_PSMs_per_gene.Genes))

  Filtered_tmp_submatrix <- tmp_submatrix[tmp_submatrix$Genes %in% Filtered_PSMs_per_gene$Genes, ]
  Filtered_tmp_submatrix <- subset(Filtered_tmp_submatrix, select=-c(Stripped.Sequence))

  Filtered_tmp_submatrix <- aggregate(Filtered_tmp_submatrix[ , 3:ncol(Filtered_tmp_submatrix)], list("Genes" = Filtered_tmp_submatrix$Genes), median, na.rm =TRUE)


  # Read Tryptic peptide file
  fasta_seq  <- read.table(file="/Users/ananth/Documents/DIANN/SwissProt_TrypticPeptides.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  seq_entry <- data.frame(table(fasta_seq$Gene))
  nondup_seq_entry <- seq_entry[seq_entry$Freq == 1, ] 

  tryptic_pept <- fasta_seq[fasta_seq$Gene %in% nondup_seq_entry$Var1, ]
  tryptic_pept <- tryptic_pept[,c("Gene","number_of_tryptic_peptides")]


  Filtered_tmp_submatrix_LFQ <- merge(x=Filtered_tmp_submatrix, y=tryptic_pept,
                                  by.x=c("Genes"), by.y=c("Gene"),
                                  all.x=FALSE, all.y=FALSE)

  # Normalise LFQs in all columns by number of tryptic peptides
  Tryptic_pept_normalised_iBAQ <- Filtered_tmp_submatrix_LFQ %>%
    mutate(
      across(c(2:(ncol(Filtered_tmp_submatrix_LFQ)-1)),
             .fns = ~./number_of_tryptic_peptides))
      # Not multiplying by a billion (no ppb scale)
      #across(c(2:(ncol(Filtered_tmp_submatrix_LFQ)-1)),
       #      .fns = ~.*1000000000))
  
  Tryptic_pept_normalised_iBAQ <- subset(Tryptic_pept_normalised_iBAQ, select = -c(number_of_tryptic_peptides) )
  colnames(Tryptic_pept_normalised_iBAQ) <- paste("iBAQ", colnames(Tryptic_pept_normalised_iBAQ),sep=".")
  colnames(Tryptic_pept_normalised_iBAQ)[1] <- "Gene"


  ##### Write normalised DIANN protein abundance matrix to export to Expression Atlas
  write.table(Tryptic_pept_normalised_iBAQ, file = paste("/Users/ananth/Documents/DIANN/",Dataset,"/SwissProt/MappedToGeneID.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE )
  
  print(paste("Finished processing dataset ....",Dataset))
}

