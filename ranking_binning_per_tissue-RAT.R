################################################
##### Files generated for OpenTargets ##########
################################################

library(ggplot2)
library(dplyr)
library(stats)
library(tidyr)
library(viridis)
library(grid)
library(RColorBrewer)
library(heatmap3)
library(stringr)
library(ggrepel)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
library(sva)

BiocManager::install("limma")
library(limma)

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/')


Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene.ID", "Gene.Symbol"), by.y=c("Gene.ID", "Gene.Symbol"), all.x=TRUE, all.y=TRUE)
}

Get_median <- function(dataset){
  dataset$Median = apply(dataset[,-c(1)], 1, median, na.rm = T)
  dataset <- dataset[,c("Gene", "Median"), drop=FALSE]
}


dataset1 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD001839/proteinGroups_ppb_final-tissue_names_PXD001839.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD003375/proteinGroups_ppb_final-tissue_names_PXD003375.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD004364/proteinGroups_ppb_final-tissue_names_PXD004364.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD006692/proteinGroups_ppb_final-tissue_names_PXD006692.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD012677/proteinGroups_ppb_final-tissue_names_PXD012677.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD013543/proteinGroups_ppb_final-tissue_names_PXD013543.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD015928/proteinGroups_ppb_final-tissue_names_PXD015928.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD016793/proteinGroups_ppb_final-tissue_names_PXD016793.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD016958/proteinGroups_ppb_final-tissue_names_PXD016958.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


Merged_input_data <- Merge_data(dataset1, dataset2)
Merged_input_data <- Merge_data(Merged_input_data, dataset3)
Merged_input_data <- Merge_data(Merged_input_data, dataset4)
Merged_input_data <- Merge_data(Merged_input_data, dataset5)
Merged_input_data <- Merge_data(Merged_input_data, dataset6)
Merged_input_data <- Merge_data(Merged_input_data, dataset7)
Merged_input_data <- Merge_data(Merged_input_data, dataset8)
Merged_input_data <- Merge_data(Merged_input_data, dataset9)

Merged_input_data[Merged_input_data == 0] <- NA

rm(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6, dataset7, dataset8, dataset9)

colnames(Merged_input_data) <- gsub("Gene.ID", "GeneID", colnames(Merged_input_data))
colnames(Merged_input_data) <- gsub("Gene.Symbol", "GeneName", colnames(Merged_input_data))
colnames(Merged_input_data) <- gsub(".*_", "", colnames(Merged_input_data), perl=TRUE)

Gene_info <- Merged_input_data[,c("GeneID", "GeneName")]

#Some of the gene entries (ex. IGHA2 has two Ensembl gene ids ENSG00000211890 & ENSG00000276173
#                          ex. IGHV2-70 has two Ensembl gene ids ENSG00000274576 & ENSG00000282453)
# because of this there are duplicate entries of such genes. These are aggregated by taking the median of them

Merged_input_data <- aggregate(Merged_input_data[,-c(1,2)], list("GeneName"= Merged_input_data$GeneName), median, na.rm =TRUE)


#2.Separate samples by tissues and extract genes
Brain_Amygdala <- Merged_input_data[grepl("GeneName|Amygdala", colnames(Merged_input_data), ignore.case = TRUE)]
Heart_LeftVentricle <-  Merged_input_data[grepl("GeneName|LeftVentricle", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_CorticalCollectingDuct <- Merged_input_data[grepl("GeneName|CorticalCollectingDuct", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_ConnectingTubule <- Merged_input_data[grepl("GeneName|ConnectingTubule", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_CorticalThickAscendingLimb <- Merged_input_data[grepl("GeneName|CorticalThickAscendingLimb", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_DistalConvolutedTubule <- Merged_input_data[grepl("GeneName|DistalConvolutedTubule", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_InnerMedullaryCollectingDuct <- Merged_input_data[grepl("GeneName|InnerMedullaryCollectingDuct", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_MedullaryThickAscendingLimb <- Merged_input_data[grepl("GeneName|MedullaryThickAscendingLimb", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_OuterMedullaryCollectingDuct <- Merged_input_data[grepl("GeneName|OuterMedullaryCollectingDuct", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_FirstSegmentOfProximalTubule <- Merged_input_data[grepl("GeneName|FirstSegmentOfProximalTubule", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_SecondSegmentOfProximalTubule <- Merged_input_data[grepl("GeneName|SecondSegmentOfProximalTubule", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney_ThirdSegmentOfProximalTubule <- Merged_input_data[grepl("GeneName|ThirdSegmentOfProximalTubule", colnames(Merged_input_data), ignore.case = TRUE)]
Liver <- Merged_input_data[grepl("GeneName|Liver", colnames(Merged_input_data), ignore.case = TRUE)]
Lung <- Merged_input_data[grepl("GeneName|Lung", colnames(Merged_input_data), ignore.case = TRUE)]
SpinalCord_CaudalSegment <- Merged_input_data[grepl("GeneName|CaudalSegment", colnames(Merged_input_data), ignore.case = TRUE)]
SpinalCord_RostralSegment <- Merged_input_data[grepl("GeneName|RostralSegment", colnames(Merged_input_data), ignore.case = TRUE)]
Tendon <- Merged_input_data[grepl("GeneName|Tendon", colnames(Merged_input_data), ignore.case = TRUE)]
Testis <- Merged_input_data[grepl("GeneName|Testis", colnames(Merged_input_data), ignore.case = TRUE)]


# Compute the median from all samples of each tissue separately
Tissue_median <- as.data.frame(cbind(GeneName = Brain_Amygdala$GeneName, Amygdala = apply(as.data.frame(Brain_Amygdala[,-c(1)]), 1, median, na.rm = T)))
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Heart_LeftVentricle$GeneName, Left_ventricle = apply(Heart_LeftVentricle[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_CorticalCollectingDuct$GeneName, Cortical_collecting_duct = apply(as.data.frame(Kidney_CorticalCollectingDuct[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_ConnectingTubule$GeneName, Connecting_tubule = apply(as.data.frame(Kidney_ConnectingTubule[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_CorticalThickAscendingLimb$GeneName, Cortical_thick_ascending_limb = apply(as.data.frame(Kidney_CorticalThickAscendingLimb[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_DistalConvolutedTubule$GeneName, Distal_convoluted_tubule = apply(as.data.frame(Kidney_DistalConvolutedTubule[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_InnerMedullaryCollectingDuct$GeneName, Inner_medullary_collecting_duct = apply(as.data.frame(Kidney_InnerMedullaryCollectingDuct[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_MedullaryThickAscendingLimb$GeneName, Medullary_thick_ascending_limb = apply(as.data.frame(Kidney_MedullaryThickAscendingLimb[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_OuterMedullaryCollectingDuct$GeneName, Outer_medullary_collecting_duct = apply(as.data.frame(Kidney_OuterMedullaryCollectingDuct[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_FirstSegmentOfProximalTubule$GeneName, First_segment_of_Proximal_tubule = apply(as.data.frame(Kidney_FirstSegmentOfProximalTubule[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_SecondSegmentOfProximalTubule$GeneName, Second_segment_of_Proximal_tubule = apply(as.data.frame(Kidney_SecondSegmentOfProximalTubule[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Kidney_ThirdSegmentOfProximalTubule$GeneName, Third_segment_of_Proximal_tubule = apply(as.data.frame(Kidney_ThirdSegmentOfProximalTubule[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Liver$GeneName, Liver = apply(Liver[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Lung$GeneName, Lung = apply(as.data.frame(Lung[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = SpinalCord_CaudalSegment$GeneName, Caudal_segment_of_Spinalcord = apply(as.data.frame(SpinalCord_CaudalSegment[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = SpinalCord_RostralSegment$GeneName, Rostral_segment_of_Spinalcord = apply(as.data.frame(SpinalCord_RostralSegment[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Tendon$GeneName, Tendon = apply(as.data.frame(Tendon[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Testis$GeneName, Testis = apply(as.data.frame(Testis[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)

write.table(Tissue_median, file = paste("Gene_distribution_in_Tissues-GeneNames-Median_intensities-RAT.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )


# Median of tissue bins
# The inout file "Binned_intensities_all_samples-RAT.txt" is obtained by reading and merging binne dexpression 
# values of tissues separately binned within a dataset.
# see script: Gene_distribution_in_organs_Median_bins-RAT.R
######
all_tissue_bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Binned_intensities_all_samples-RAT.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

all_tissue_bins_long <- gather(all_tissue_bins, Datasets, Bins, colnames(all_tissue_bins[2]):colnames(all_tissue_bins[ncol(all_tissue_bins)]), factor_key=TRUE)
all_tissue_bins_long$Tissues <- gsub(".*\\.", "", all_tissue_bins_long$Datasets)
all_tissue_bins_long$Tissues <- gsub("CaudalSegment.*", "CaudalSegment_of_SpinalCord", all_tissue_bins_long$Tissues)
all_tissue_bins_long$Tissues <- gsub("RostralSegment.*", "RostralSegment_of_SpinalCord", all_tissue_bins_long$Tissues)

all_tissue_bins_long_aggregate <- aggregate(all_tissue_bins_long[,-c(2)], list("Gene" = all_tissue_bins_long$Gene, "Tissues" = all_tissue_bins_long$Tissues), median, na.rm =TRUE)

all_tissue_bins_long_aggregate <- all_tissue_bins_long_aggregate[,-c(3,5)]


All_tissues_median_bins<- spread(all_tissue_bins_long_aggregate, Tissues, Bins)

write.table(All_tissues_median_bins, file = "Gene_distribution_in_Tissues-GeneNames_bin_values-RAT.txt", sep = "\t", row.names = FALSE, quote = FALSE )


##############################
#### To double check the "convert_to_ranked_data" method as the output only has one instance of bin1
##############################



#foo1 <- data.frame(GeneName=c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG"), A=c(10,NA,3,4,10,20,30), B=c(10,20,30,40,50,60,70), C=c(70,12,23,33,44,55,65))
#gene_names <- foo1[,c(1)]

gene_names <- Tissue_median[,c(1)] 


convert_to_scaled_data <- function(data, bins) {
  apply(data, 1, function(sample) {
    ceiling((sample / max(sample, na.rm = T)) * bins)
  })
}

convert_to_ranked_data <- function(data) {
  apply(data, 2, function(sample_data){
    unique_vals = unique(sample_data)
    #unique_vals = unique_vals[order(unique_vals)]
    unique_vals = unique_vals[order(unique_vals, na.last = NA)]
    #sapply(sample_data, function(s){which(s == unique_vals)}) - 1
    sapply(sample_data, function(s){ if(!is.na(s)){which(s == unique_vals)} else {s=NA} }) - 1
  })
}

Binning <- function(data_input, gene_names) {
  # Make a data frame of sample + protein/peptide abundances. 
  data = data.frame(data_input)
  
  # Rank data.
  data_ranked = as.data.frame(t(convert_to_ranked_data(data)))
  colnames(data_ranked) <- gene_names
  
  # Different size of bins to try.
  binning_size_approaches = list(
    bins_2 = 1,
    bins_5 = 4,
    bins_10 = 9,
    bins_100 = 99,
    bins_1000 = 999,
    bins_10000 = 9999
  )
  
  # Transform data using different bin sizes.
  results = lapply(binning_size_approaches, function(binning_size_approach) {
    # 'Bin' data by scaling the ranked data.
    binned_data = as.data.frame(t(convert_to_scaled_data(data_ranked, binning_size_approach)))
  })
}

# If a dataset has different tissues or tissue regions, do not bin all together,
# these tissues/regions must be separated and binned separately.

colnames(Tissue_median)

binning_input <- Tissue_median[,-c(1)]

#backup before changing NAs to 0
#binning_input_with_NA <- binning_input


# non-detected genes (NA or 0 values) are assigned 0 values here to help in binning.
#binning_input[is.na(binning_input)] <- 0


#Binning
#foo_bins <- Binning(foo1[,-c(1)], gene_names)
Dataset_all_bins <- Binning(binning_input, gene_names)


# Consider only data grouped into 5 bins
Get_Bin_5 <- function(binned_data_list){
  
  Bins_5 <- as.data.frame(t(binned_data_list$bins_5))
  # rename bins 0 to 4 --> 1 to 5
  # lowest intensities are in bin1 and highest in bin5
  Bins_5[Bins_5 == 4] <- 5
  Bins_5[Bins_5 == 3] <- 4
  Bins_5[Bins_5 == 2] <- 3
  Bins_5[Bins_5 == 1] <- 2
  Bins_5[Bins_5 == 0] <- 1
  
  Bins_5 <- tibble::rownames_to_column(Bins_5, "Gene")
}

Dataset_5_bins <- Get_Bin_5(Dataset_all_bins)


# Now reassign NA values back to those genes which were undetected or 0 in ALL samples 
# NOTE: Only those genes are reassigned back to NA that have missing values in
# all samples (ie., if at least one sample had a value then it is left as 0)
# this particular reassigning is done after binning (later in the code).


#Dataset_5_bins[which(apply(binning_input_with_NA, 1, function(x) all(is.na(x)) )),-c(1)] <- NA


write.table(Dataset_5_bins, file = "Gene_distribution_in_Tissues-GeneNames_bin_values-TORECHECK_RAT.txt", sep = "\t", row.names = FALSE, quote = FALSE )

