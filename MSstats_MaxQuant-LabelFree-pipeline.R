# MSstats tutorial for post-processing MaxQuant outputs - LabelFree
# to obtain protein quantification (abundance)

# http://msstats.org/wp-content/uploads/2019/11/MSstats_v3.18.1_manual.pdf


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats")


library('MSstats', warn.conflicts = F, quietly = T, verbose = F)

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD005819_33threads_yoda/')
# Tutorial example datasets here: ftp://massive.ucsd.edu/RMSV000000249/2019-06-03_mnchoi_8befbb21

# 1. Read proteinGroups to get proteinID information.
# proteinGroups  <- read.table( "/Users/ananth/Documents/MSstats/Example_datasets/Choi2017_DDA_MaxQuant_proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
proteinGroups  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


# 2. Read evidence.txt file
# evidence_file  <- read.table( "/Users/ananth/Documents/MSstats/Example_datasets/Choi2017_DDA_MaxQuant_evidence.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
evidence_file  <- read.table( "evidence.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


# 3. Read annotation file
# Note In this file the RAW file names should not have .raw extension!
# annot <- read.table( "/Users/ananth/Documents/MSstats/Example_datasets/Choi2017_DDA_MaxQuant_annotation.csv" , quote = "\"", header = TRUE, sep = ",", stringsAsFactors = FALSE, comment.char = "#")
annot <- read.table( "MSstat_annot.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


# 4. Convert MaxQuant output file into format required for MSstats
quant <- MaxQtoMSstatsFormat(evidence=evidence_file, annotation=annot, proteinGroups=proteinGroups, removeProtein_with1Peptide=TRUE)


# 5. Normalisation step
maxquant.processed <- dataProcess(quant, normalization='equalizeMedians',
                                 summaryMethod="TMP", 
                                 cutoffCensored="minFeature",
                                 MBimpute = TRUE,
                                 censoredInt="NA", ## !! important for MaxQuant MBimpute=TRUE,
                                 maxQuantileforCensored=0.999)


# names(maxquant.processed)
# head(maxquant.processed$ProcessedData)
# head(maxquant.processed$RunlevelData)
# levels(maxquant.processed$ProcessedData$GROUP_ORIGINAL)

# 6. Protein quantification in each group
# https://rdrr.io/bioc/MSstats/man/quantification.html
Protein_quantify <- quantification(maxquant.processed, type = "Group", format = "matrix")


write.table(Protein_quantify, "Protein_abundance_MSstats.txt", sep = "\t", row.names = FALSE, quote = FALSE )

