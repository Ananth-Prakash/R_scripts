#Rebinned using Excel. on median intensities files
# Aggregate protein groups

#Rebinned to 5 bins
Excel_rebinned_bins5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/DDA_Manuscript_Supplementary_Files/Rebinned_5_bins.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
Excel_rebinned_bins5 <- Excel_rebinned_bins5[,-c(1)]

Aggregate_Excel_rebinned_bins5 <-  aggregate(Excel_rebinned_bins5[ , 2:(ncol(Excel_rebinned_bins5))], list("Gene Symbol" = Excel_rebinned_bins5$Gene.Symbol, "Ensembl.ID" = Excel_rebinned_bins5$EnsemblID), median, na.rm =TRUE)
Aggregate_Excel_rebinned_bins5 <- Aggregate_Excel_rebinned_bins5[,-c(3)]

Gene_info <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/DDA_Manuscript_Supplementary_Files/IDs.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


Aggregate_Excel_rebinned_bins5 <- merge(x=Gene_info, y=Aggregate_Excel_rebinned_bins5,
                                 by.x=c("GeneName"), by.y=c("Gene Symbol"))

Aggregate_Excel_rebinned_bins5 <- Aggregate_Excel_rebinned_bins5[,-c(4)]
Aggregate_Excel_rebinned_bins5 <-  aggregate(Aggregate_Excel_rebinned_bins5[ , 4:(ncol(Aggregate_Excel_rebinned_bins5))], list("Gene.Name" = Aggregate_Excel_rebinned_bins5$GeneName, "Gene.ID" = Aggregate_Excel_rebinned_bins5$GeneID, "UniProt.ID" = Aggregate_Excel_rebinned_bins5$UniProt_ID), median, na.rm =TRUE)


write.table(Aggregate_Excel_rebinned_bins5, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/DDA_Manuscript_Supplementary_Files/Aggregate_Excel_rebinned_bins5.txt", sep = "\t", row.names = FALSE, quote = FALSE )
  

#Rebinned to 3 bins
Excel_rebinned_bins3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/DDA_Manuscript_Supplementary_Files/Rebinned_3_bins.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
Excel_rebinned_bins3 <- Excel_rebinned_bins3[,-c(1)]

Aggregate_Excel_rebinned_bins3 <-  aggregate(Excel_rebinned_bins3[ , 2:(ncol(Excel_rebinned_bins3))], list("Gene Symbol" = Excel_rebinned_bins3$Gene.Symbol, "Ensembl.ID" = Excel_rebinned_bins3$EnsemblID), median, na.rm =TRUE)
Aggregate_Excel_rebinned_bins3 <- Aggregate_Excel_rebinned_bins3[,-c(3)]

Aggregate_Excel_rebinned_bins3 <- merge(x=Gene_info, y=Aggregate_Excel_rebinned_bins3,
                                        by.x=c("GeneName"), by.y=c("Gene Symbol"))

Aggregate_Excel_rebinned_bins3 <- Aggregate_Excel_rebinned_bins3[,-c(4)]
Aggregate_Excel_rebinned_bins3 <-  aggregate(Aggregate_Excel_rebinned_bins3[ , 4:(ncol(Aggregate_Excel_rebinned_bins3))], list("Gene.Name" = Aggregate_Excel_rebinned_bins3$GeneName, "Gene.ID" = Aggregate_Excel_rebinned_bins3$GeneID, "UniProt.ID" = Aggregate_Excel_rebinned_bins3$UniProt_ID), median, na.rm =TRUE)



write.table(Aggregate_Excel_rebinned_bins3, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/DDA_Manuscript_Supplementary_Files/Aggregate_Excel_rebinned_bins3.txt", sep = "\t", row.names = FALSE, quote = FALSE )

