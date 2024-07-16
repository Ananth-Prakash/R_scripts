# To get the median bin values of organs as Summary after binning tissues separately
library(RColorBrewer)
library(heatmap3)
library(pheatmap)
library(Hmisc)
setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/')

##### Secretome data
CRC_Secretome <- read.table("CRC_Secretome/proteinGroups_ppb_final_binning.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# Keep the minimum number of genes that are common in all datsets
#Merge_data <- function(dataset1, dataset2){
#  merged <- merge(x=dataset1, y=dataset2,
#                  by.x=c("Gene"), by.y=c("Gene"), all.x=TRUE, all.y=TRUE)
#}

#merged_data <- Merge_data(dataset1, dataset2)


#binned_samples <- merged_data
#binned_samples <- CRC_Secretome

cor_results1 <- rcorr(as.matrix(CRC_Secretome[,-c(1,2)]), type = c("spearman"))
correl_matrix <- cor_results1$r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correl_summary <- flattenCorrMatrix(cor_results1$r, cor_results1$P) 



Pairwise_only_bloodderived <- with(correl_summary, correl_summary[ grepl("Blood.derived", row) & grepl("Blood.derived", column),])
Pairwise_only_bloodderived$R2 <- Pairwise_only_bloodderived$cor*Pairwise_only_bloodderived$cor

Pairwise_only_cellculture <- with(correl_summary, correl_summary[ grepl("Cellculture.derived", row) & grepl("Cellculture.derived", column),])
Pairwise_only_cellculture$R2 <- Pairwise_only_cellculture$cor*Pairwise_only_cellculture$cor

Pairwise_only_InterstitialFluid_Normal <- with(correl_summary, correl_summary[ grepl("InterstitialFluid.Normal", row) & grepl("InterstitialFluid.Normal", column),])
Pairwise_only_InterstitialFluid_Normal$R2 <- Pairwise_only_InterstitialFluid_Normal$cor*Pairwise_only_InterstitialFluid_Normal$cor

Pairwise_only_InterstitialFluid_ColonCancer <- with(correl_summary, correl_summary[ grepl("InterstitialFluid.ColorectalCancer", row) & grepl("InterstitialFluid.ColorectalCancer", column),])
Pairwise_only_InterstitialFluid_ColonCancer$R2 <- Pairwise_only_InterstitialFluid_ColonCancer$cor*Pairwise_only_InterstitialFluid_ColonCancer$cor

Pairwise_only_ECvesicles_and_exosomes_Normal <- with(correl_summary, correl_summary[ grepl("ECvesicles_and_exosomes.Normal", row) & grepl("ECvesicles_and_exosomes.Normal", column),])
Pairwise_only_ECvesicles_and_exosomes_Normal$R2 <- Pairwise_only_ECvesicles_and_exosomes_Normal$cor*Pairwise_only_ECvesicles_and_exosomes_Normal$cor

Pairwise_only_ECvesicles_and_exosomes_ColonCancer <- with(correl_summary, correl_summary[ grepl("ECvesicles_and_exosomes.ColorectalCancer", row) & grepl("ECvesicles_and_exosomes.ColorectalCancer", column),])
Pairwise_only_ECvesicles_and_exosomes_ColonCancer$R2 <- Pairwise_only_ECvesicles_and_exosomes_ColonCancer$cor*Pairwise_only_ECvesicles_and_exosomes_ColonCancer$cor

Significant <- correl_summary[correl_summary$p < 0.05,]

#Organs <- data.frame(ID= factor(gsub(".*Sample[0-9]+\\.", "", colnames(binned_samples)[-c(1,2)], perl=TRUE)))

annotation <- data.frame(Samples = gsub(".*Sample[0-9]+\\.", "", colnames(correl_matrix), perl=TRUE))
annotation$Samples <- gsub("\\.ColorectalCancer","", annotation$Samples)
annotation$Samples <- gsub("\\.Normal","", annotation$Samples)
annotation$Samples <- gsub("Blood\\.derived","Blood-derived", annotation$Samples)
annotation$Samples <- gsub("Cellculture\\.derived","Cell culture-derived", annotation$Samples)
annotation$Samples <- gsub("ECvesicles_and_exosomes","EC vesicles and exosomes", annotation$Samples)
annotation$Samples <- gsub("InterstitialFluid","Interstitial fluid", annotation$Samples)

rownames(annotation) <- colnames(correl_matrix)
annotation$Condition <- gsub(".*\\.", "", colnames(correl_matrix), perl=TRUE)
annotation$Condition <- gsub("ColorectalCancer","Colorectal cancer", annotation$Condition)
annotation$Datasets <- gsub("\\..*", "", rownames(annotation), perl=TRUE)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation$Samples))))
annoCol <- newCols(length(unique(annotation$Samples)))
names(annoCol) <- unique(annotation$Samples)
annoCol <- list(category = annoCol)

ann_colors = list(
  Datasets = c(JPST000867="#E69F00", PXD005693="#009E73", PXD005709="#F0E442",
               PXD010458="#0072B2", PXD020454="#FF595E", PXD031556="#CC79A7",
               PXD032899="#000000"),
  Samples = c("Blood-derived"="#FF595E","Cell culture-derived"="#FFCA3A",
              "EC vesicles and exosomes"="#69BE28", "Interstitial fluid"="#1982C4"))

#PXD020454="#D55E00"
pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
         annotation_col = annotation,
         border_color = NA,
         #annotation_colors = Sample_colours[1],
         #annotation_colors = annoCol,
         annotation_colors = ann_colors,
         annotation = annotation, 
         fontsize = 6)

#Clear previous data
rm(list=ls()) 

######## Tumor data
CRC_Tumor     <- read.table("CRC_Tumour/proteinGroups_ppb_final_binning.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(CRC_Tumor) <- gsub("Normal","Mucosa", colnames(CRC_Tumor))

cor_results1 <- rcorr(as.matrix(CRC_Tumor[,-c(1,2)]), type = c("spearman"))
correl_matrix <- cor_results1$r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correl_summary <- flattenCorrMatrix(cor_results1$r, cor_results1$P) 



Pairwise_only_Tumor <- with(correl_summary, correl_summary[ grepl("Tumor", row) & grepl("Tumor", column),])
Pairwise_only_Tumor$R2 <- Pairwise_only_Tumor$cor*Pairwise_only_Tumor$cor

Pairwise_only_Mucosa <- with(correl_summary, correl_summary[ grepl("Mucosa", row) & grepl("Mucosa", column),])
Pairwise_only_Mucosa$R2 <- Pairwise_only_Mucosa$cor*Pairwise_only_Mucosa$cor

Pairwise_only_Adenoma <- with(correl_summary, correl_summary[ grepl("Adenoma", row) & grepl("Adenoma", column),])
Pairwise_only_Adenoma$R2 <- Pairwise_only_Adenoma$cor*Pairwise_only_Adenoma$cor

Significant <- correl_summary[correl_summary$p < 0.05,]

#Organs <- data.frame(ID= factor(gsub(".*Sample[0-9]+\\.", "", colnames(binned_samples)[-c(1,2)], perl=TRUE)))

annotation <- data.frame(Samples = gsub(".*Sample[0-9]+\\.", "", colnames(correl_matrix), perl=TRUE))
annotation$Samples <- gsub("\\..*","", annotation$Samples)
annotation$Samples <- gsub("Colonmucosa","Colon mucosa", annotation$Samples)

rownames(annotation) <- colnames(correl_matrix)
annotation$Condition <- gsub(".*\\.", "", colnames(correl_matrix), perl=TRUE)
annotation$Datasets <- gsub("\\..*", "", rownames(annotation), perl=TRUE)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation$Samples))))
annoCol <- newCols(length(unique(annotation$Samples)))
names(annoCol) <- unique(annotation$Samples)
annoCol <- list(category = annoCol)

annotation <- subset(annotation, select=-c(Samples))

ann_colors = list(
  Datasets = c(CPTAC="#E69F00", PXD001676="#009E73", PXD002137="#F0E442",
               PXD014511="#0072B2", PXD019504="#FF595E"),
  #Samples = c("Colon"="#FF595E","Colon mucosa"="#FFCA3A",
  #            "Colorectum"="#69BE28"),
  Condition = c("Adenoma"="#17BEBB","Mucosa"="#A73030FF", "Tumor"="#000000"))

#PXD020454="#D55E00"
pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
         annotation_col = annotation,
         border_color = NA,
         #annotation_colors = Sample_colours[1],
         #annotation_colors = annoCol,
         annotation_colors = ann_colors,
         annotation = annotation, 
         fontsize = 6)

