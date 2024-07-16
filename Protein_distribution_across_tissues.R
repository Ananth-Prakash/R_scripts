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

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/')


Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene.ID", "Gene.Name"), by.y=c("Gene.ID", "Gene.Name"), all.x=TRUE, all.y=TRUE)
}

Get_median <- function(dataset){
  dataset$Median = apply(dataset[,-c(1)], 1, median, na.rm = T)
  dataset <- dataset[,c("Gene", "Median"), drop=FALSE]
}


dataset1 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000547/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000548/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001325/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001608_30threads_yoda/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD002029/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004143/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004332/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD005819_33threads_yoda/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006233/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset10 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006675/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset11 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008934/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset12 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010154_Ananth/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset13 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010271/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset14 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012131/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset15 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012755/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset16 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD015079/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset17 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD020187/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset18 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/ACT_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset19 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset20 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Banner_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset21 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset22 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_Precuneus/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset23 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Mayo_TemporalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset24 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset25 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset26 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012431/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset27 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008722/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")



colnames(dataset12)[2] <- "Gene.Name"
colnames(dataset14)[2] <- "Gene.Name"
colnames(dataset16)[2] <- "Gene.Name"
colnames(dataset17)[2] <- "Gene.Name"
colnames(dataset18)[2] <- "Gene.Name"
colnames(dataset19)[2] <- "Gene.Name"
colnames(dataset20)[2] <- "Gene.Name"
colnames(dataset21)[2] <- "Gene.Name"
colnames(dataset22)[2] <- "Gene.Name"
colnames(dataset23)[2] <- "Gene.Name"
colnames(dataset25)[2] <- "Gene.Name"


Merged_input_data <- Merge_data(dataset1, dataset2)
Merged_input_data <- Merge_data(Merged_input_data, dataset3)
Merged_input_data <- Merge_data(Merged_input_data, dataset4)
Merged_input_data <- Merge_data(Merged_input_data, dataset5)
Merged_input_data <- Merge_data(Merged_input_data, dataset6)
Merged_input_data <- Merge_data(Merged_input_data, dataset7)
Merged_input_data <- Merge_data(Merged_input_data, dataset8)
Merged_input_data <- Merge_data(Merged_input_data, dataset9)
Merged_input_data <- Merge_data(Merged_input_data, dataset10)
Merged_input_data <- Merge_data(Merged_input_data, dataset11)
Merged_input_data <- Merge_data(Merged_input_data, dataset12)
Merged_input_data <- Merge_data(Merged_input_data, dataset13)
Merged_input_data <- Merge_data(Merged_input_data, dataset14)
Merged_input_data <- Merge_data(Merged_input_data, dataset15)
Merged_input_data <- Merge_data(Merged_input_data, dataset16)
Merged_input_data <- Merge_data(Merged_input_data, dataset17)
Merged_input_data <- Merge_data(Merged_input_data, dataset18)
Merged_input_data <- Merge_data(Merged_input_data, dataset19)
Merged_input_data <- Merge_data(Merged_input_data, dataset20)
Merged_input_data <- Merge_data(Merged_input_data, dataset21)
Merged_input_data <- Merge_data(Merged_input_data, dataset22)
Merged_input_data <- Merge_data(Merged_input_data, dataset23)
Merged_input_data <- Merge_data(Merged_input_data, dataset24)
Merged_input_data <- Merge_data(Merged_input_data, dataset25)
#Merged_input_data <- Merge_data(Merged_input_data, dataset26)
#Merged_input_data <- Merge_data(Merged_input_data, dataset27)

Merged_input_data[Merged_input_data == 0] <- NA

rm(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6, dataset7, dataset8, dataset9, dataset10, 
   dataset11, dataset12, dataset13, dataset14, dataset15, dataset16, dataset17, dataset18, dataset19, dataset20, dataset21, dataset22, dataset23, dataset24, dataset25)

colnames(Merged_input_data) <- gsub("Gene.ID", "GeneID", colnames(Merged_input_data))
colnames(Merged_input_data) <- gsub("Gene.Name", "GeneName", colnames(Merged_input_data))
colnames(Merged_input_data) <- gsub(".*_", "", colnames(Merged_input_data), perl=TRUE)
#colnames(Merged_input_data) <- gsub(".*\\.", "", colnames(Merged_input_data), perl=TRUE)
colnames(Merged_input_data) <- gsub("breast", "Breast", colnames(Merged_input_data), perl=TRUE)
colnames(Merged_input_data) <- gsub("UmblicalArtery", "UmbilicalArtery", colnames(Merged_input_data), perl=TRUE )

Gene_info <- Merged_input_data[,c("GeneID", "GeneName")]
#Gene_info <- Gene_info[!duplicated(Gene_info$GeneID), ]

#Some of the gene entries (ex. IGHA2 has two Ensembl gene ids ENSG00000211890 & ENSG00000276173
#                          ex. IGHV2-70 has two Ensembl gene ids ENSG00000274576 & ENSG00000282453)
# because of this there are duplicate entries of such genes. These are aggregated by taking the median of them

foo <- Merged_input_data
gene_maps <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_UniProt_IDs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
foo <- merge(x=gene_maps, y=Merged_input_data,
             by.x=c("Gene.ID","Gene.Name"), by.y=c("GeneID","GeneName"),
             all.x=FALSE, all.y=TRUE)
write.table(foo, file = paste("Ppb_intensities_all_samples.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

#all_bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Binned_intensities_all_samples.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#all_bins <- merge(x=gene_maps, y=all_bins,
#                  by.x=c("Gene.Name"), by.y=c("Gene"),
#                  a..x=FALSE, all.y=TRUE)
#write.table(all_bins, file = paste("binned_intensities_all_samples.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )


Merged_input_data <- aggregate(Merged_input_data[,-c(1,2)], list("GeneName"= Merged_input_data$GeneName), median, na.rm =TRUE)


#Separate samples by organs and extract genes
AdiposeTissue <- Merged_input_data[grepl("GeneName|AdiposeTissue", colnames(Merged_input_data), ignore.case = TRUE)]
AdrenalGland <- Merged_input_data[grepl("GeneName|AdrenalGland", colnames(Merged_input_data), ignore.case = TRUE)]
BoneMarrow <- Merged_input_data[grepl("GeneName|BoneMarrow", colnames(Merged_input_data), ignore.case = TRUE)]
Brain <- Merged_input_data[grepl("GeneName|Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PituitaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex", colnames(Merged_input_data), ignore.case = TRUE)]
#Breast <- Merged_input_data[grepl("GeneName|Breast", colnames(Merged_input_data), ignore.case = TRUE)]
Colon <- Merged_input_data[grepl("GeneName|Colon", colnames(Merged_input_data), ignore.case = TRUE)]
Duodenum <- Merged_input_data[grepl("GeneName|Duodenum", colnames(Merged_input_data), ignore.case = TRUE)]
Esophagus <- Merged_input_data[grepl("GeneName|Esophagus", colnames(Merged_input_data), ignore.case = TRUE)]
FallopianTubeOviduct <- Merged_input_data[grepl("GeneName|FallopianTubeOviduct", colnames(Merged_input_data), ignore.case = TRUE)]
GallBladder <- Merged_input_data[grepl("GeneName|GallBladder", colnames(Merged_input_data), ignore.case = TRUE)]
Heart <-  Merged_input_data[grepl("GeneName|Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney <- Merged_input_data[grepl("GeneName|Kidney", colnames(Merged_input_data), ignore.case = TRUE)]
Liver <- Merged_input_data[grepl("GeneName|Liver", colnames(Merged_input_data), ignore.case = TRUE)]
Lung <- Merged_input_data[grepl("GeneName|Lung", colnames(Merged_input_data), ignore.case = TRUE)]
LymphNode <- Merged_input_data[grepl("GeneName|LymphNode", colnames(Merged_input_data), ignore.case = TRUE)]
Ovary <- Merged_input_data[grepl("GeneName|Ovary", colnames(Merged_input_data), ignore.case = TRUE)]
Pancreas <- Merged_input_data[grepl("GeneName|Pancreas", colnames(Merged_input_data), ignore.case = TRUE)]
Placenta <- Merged_input_data[grepl("GeneName|Placenta", colnames(Merged_input_data), ignore.case = TRUE)]
Prostate <- Merged_input_data[grepl("GeneName|Prostate", colnames(Merged_input_data), ignore.case = TRUE)]
Rectum <- Merged_input_data[grepl("GeneName|Rectum", colnames(Merged_input_data), ignore.case = TRUE)]
SalivaryGland <- Merged_input_data[grepl("GeneName|SalivaryGland", colnames(Merged_input_data), ignore.case = TRUE)]
SmallIntestine <- Merged_input_data[grepl("GeneName|SmallIntestine", colnames(Merged_input_data), ignore.case = TRUE)]
SmoothMuscle <- Merged_input_data[grepl("GeneName|SmoothMuscle", colnames(Merged_input_data), ignore.case = TRUE)]
Spleen <- Merged_input_data[grepl("GeneName|Spleen", colnames(Merged_input_data), ignore.case = TRUE)]
Stomach <- Merged_input_data[grepl("GeneName|Stomach", colnames(Merged_input_data), ignore.case = TRUE)]
Testis <- Merged_input_data[grepl("GeneName|Testis", colnames(Merged_input_data), ignore.case = TRUE)]
Thyroid <- Merged_input_data[grepl("GeneName|Thyroid", colnames(Merged_input_data), ignore.case = TRUE)]
Tonsil <- Merged_input_data[grepl("GeneName|Tonsil", colnames(Merged_input_data), ignore.case = TRUE)]
UmbilicalArtery <- Merged_input_data[grepl("GeneName|UmbilicalArtery", colnames(Merged_input_data), ignore.case = TRUE)]
UrinaryBladder <- Merged_input_data[grepl("GeneName|UrinaryBladder", colnames(Merged_input_data), ignore.case = TRUE)]
UterineEndometrium <- Merged_input_data[grepl("GeneName|UterineEndometrium", colnames(Merged_input_data), ignore.case = TRUE)]
VermiformAppendix <- Merged_input_data[grepl("GeneName|VermiformAppendix", colnames(Merged_input_data), ignore.case = TRUE)]

# Get sample sizes for each organ (i.e., number of MS runs representing each organ)
sample_sizes <- as.data.frame(ncol(AdiposeTissue)-1)
sample_sizes <- cbind(sample_sizes, ncol(AdrenalGland)-1)
sample_sizes <- cbind(sample_sizes, ncol(BoneMarrow)-1)
sample_sizes <- cbind(sample_sizes, ncol(Brain)-1)
sample_sizes <- cbind(sample_sizes, ncol(Breast)-1)
sample_sizes <- cbind(sample_sizes, ncol(Colon)-1)
sample_sizes <- cbind(sample_sizes, ncol(Duodenum)-1)
sample_sizes <- cbind(sample_sizes, ncol(Esophagus)-1)
sample_sizes <- cbind(sample_sizes, ncol(FallopianTubeOviduct)-1)
sample_sizes <- cbind(sample_sizes, ncol(GallBladder)-1)
sample_sizes <- cbind(sample_sizes, ncol(Heart)-1)
sample_sizes <- cbind(sample_sizes, ncol(Kidney)-1)
sample_sizes <- cbind(sample_sizes, ncol(Liver)-1)
sample_sizes <- cbind(sample_sizes, ncol(Lung)-1)
sample_sizes <- cbind(sample_sizes, ncol(LymphNode)-1)
sample_sizes <- cbind(sample_sizes, ncol(Ovary)-1)
sample_sizes <- cbind(sample_sizes, ncol(Pancreas)-1)
sample_sizes <- cbind(sample_sizes, ncol(Placenta)-1)
sample_sizes <- cbind(sample_sizes, ncol(Prostate)-1)
sample_sizes <- cbind(sample_sizes, ncol(Rectum)-1)
sample_sizes <- cbind(sample_sizes, ncol(SalivaryGland)-1)
sample_sizes <- cbind(sample_sizes, ncol(SmallIntestine)-1)
sample_sizes <- cbind(sample_sizes, ncol(SmoothMuscle)-1)
sample_sizes <- cbind(sample_sizes, ncol(Spleen)-1)
sample_sizes <- cbind(sample_sizes, ncol(Stomach)-1)
sample_sizes <- cbind(sample_sizes, ncol(Testis)-1)
sample_sizes <- cbind(sample_sizes, ncol(Thyroid)-1)
sample_sizes <- cbind(sample_sizes, ncol(Tonsil)-1)
sample_sizes <- cbind(sample_sizes, ncol(UmbilicalArtery)-1)
sample_sizes <- cbind(sample_sizes, ncol(UrinaryBladder)-1)
sample_sizes <- cbind(sample_sizes, ncol(UterineEndometrium)-1)
sample_sizes <- cbind(sample_sizes, ncol(VermiformAppendix)-1)

colnames(sample_sizes) <- gsub(".*\\(", "", colnames(sample_sizes))
colnames(sample_sizes) <- gsub("\\).*", "", colnames(sample_sizes))
sample_sizes <- as.data.frame(t(sample_sizes))
sample_sizes <- tibble::rownames_to_column(sample_sizes, "Organs")
colnames(sample_sizes) <- c("Organs", "Number_of_samples")

# Data to plot distribution of iBAQ values across organs
AdiposeTissue_ppb_long <- gather(AdiposeTissue, Sample, ppb, colnames(AdiposeTissue)[2]:colnames(AdiposeTissue)[ncol(AdiposeTissue)], factor_key=TRUE)
AdrenalGland_ppb_long <- gather(AdrenalGland, Sample, ppb, colnames(AdrenalGland)[2]:colnames(AdrenalGland)[ncol(AdrenalGland)], factor_key=TRUE)
BoneMarrow_ppb_long <- gather(BoneMarrow, Sample, ppb, colnames(BoneMarrow)[2]:colnames(BoneMarrow)[ncol(BoneMarrow)], factor_key=TRUE)
Brain_ppb_long <- gather(Brain, Sample, ppb, colnames(Brain)[2]:colnames(Brain)[ncol(Brain)], factor_key=TRUE)
#Breast_ppb_long <- gather(Breast, Sample, ppb, colnames(Breast)[2]:colnames(Breast)[ncol(Breast)], factor_key=TRUE)
Colon_ppb_long <- gather(Colon, Sample, ppb, colnames(Colon)[2]:colnames(Colon)[ncol(Colon)], factor_key=TRUE)
Duodenum_ppb_long <- gather(Duodenum, Sample, ppb, colnames(Duodenum)[2]:colnames(Duodenum)[ncol(Duodenum)], factor_key=TRUE)
Esophagus_ppb_long <- gather(Esophagus, Sample, ppb, colnames(Esophagus)[2]:colnames(Esophagus)[ncol(Esophagus)], factor_key=TRUE)
FallopianTubeOviduct_ppb_long <- gather(FallopianTubeOviduct, Sample, ppb, colnames(FallopianTubeOviduct)[2]:colnames(FallopianTubeOviduct)[ncol(FallopianTubeOviduct)], factor_key=TRUE)
GallBladder_ppb_long <- gather(GallBladder, Sample, ppb, colnames(GallBladder)[2]:colnames(GallBladder)[ncol(GallBladder)], factor_key=TRUE)
Heart_ppb_long <- gather(Heart, Sample, ppb, colnames(Heart)[2]:colnames(Heart)[ncol(Heart)], factor_key=TRUE)
Kidney_ppb_long <- gather(Kidney, Sample, ppb, colnames(Kidney)[2]:colnames(Kidney)[ncol(Kidney)], factor_key=TRUE)
Liver_ppb_long <- gather(Liver, Sample, ppb, colnames(Liver)[2]:colnames(Liver)[ncol(Liver)], factor_key=TRUE)
Lung_ppb_long <- gather(Lung, Sample, ppb, colnames(Lung)[2]:colnames(Lung)[ncol(Lung)], factor_key=TRUE)
LymphNode_ppb_long <- gather(LymphNode, Sample, ppb, colnames(LymphNode)[2]:colnames(LymphNode)[ncol(LymphNode)], factor_key=TRUE)
Ovary_ppb_long <- gather(Ovary, Sample, ppb, colnames(Ovary)[2]:colnames(Ovary)[ncol(Ovary)], factor_key=TRUE)
Pancreas_ppb_long <- gather(Pancreas, Sample, ppb, colnames(Pancreas)[2]:colnames(Pancreas)[ncol(Pancreas)], factor_key=TRUE)
Placenta_ppb_long <- gather(Placenta, Sample, ppb, colnames(Placenta)[2]:colnames(Placenta)[ncol(Placenta)], factor_key=TRUE)
Prostate_ppb_long <- gather(Prostate, Sample, ppb, colnames(Prostate)[2]:colnames(Prostate)[ncol(Prostate)], factor_key=TRUE)
Rectum_ppb_long <- gather(Rectum, Sample, ppb, colnames(Rectum)[2]:colnames(Rectum)[ncol(Rectum)], factor_key=TRUE)
SalivaryGland_ppb_long <- gather(SalivaryGland, Sample, ppb, colnames(SalivaryGland)[2]:colnames(SalivaryGland)[ncol(SalivaryGland)], factor_key=TRUE)
SmallIntestine_ppb_long <- gather(SmallIntestine, Sample, ppb, colnames(SmallIntestine)[2]:colnames(SmallIntestine)[ncol(SmallIntestine)], factor_key=TRUE)
SmoothMuscle_ppb_long <- gather(SmoothMuscle, Sample, ppb, colnames(SmoothMuscle)[2]:colnames(SmoothMuscle)[ncol(SmoothMuscle)], factor_key=TRUE)
Spleen_ppb_long <- gather(Spleen, Sample, ppb, colnames(Spleen)[2]:colnames(Spleen)[ncol(Spleen)], factor_key=TRUE)
Stomach_ppb_long <- gather(Stomach, Sample, ppb, colnames(Stomach)[2]:colnames(Stomach)[ncol(Stomach)], factor_key=TRUE)
Testis_ppb_long <- gather(Testis, Sample, ppb, colnames(Testis)[2]:colnames(Testis)[ncol(Testis)], factor_key=TRUE)
Thyroid_ppb_long <- gather(Thyroid, Sample, ppb, colnames(Thyroid)[2]:colnames(Thyroid)[ncol(Thyroid)], factor_key=TRUE)
Tonsil_ppb_long <- gather(Tonsil, Sample, ppb, colnames(Tonsil)[2]:colnames(Tonsil)[ncol(Tonsil)], factor_key=TRUE)
UmbilicalArtery_ppb_long <- gather(UmbilicalArtery, Sample, ppb, colnames(UmbilicalArtery)[2]:colnames(UmbilicalArtery)[ncol(UmbilicalArtery)], factor_key=TRUE)
UrinaryBladder_ppb_long <- gather(UrinaryBladder, Sample, ppb, colnames(UrinaryBladder)[2]:colnames(UrinaryBladder)[ncol(UrinaryBladder)], factor_key=TRUE)
UterineEndometrium_ppb_long <- gather(UterineEndometrium, Sample, ppb, colnames(UterineEndometrium)[2]:colnames(UterineEndometrium)[ncol(UterineEndometrium)], factor_key=TRUE)
VermiformAppendix_ppb_long <- gather(VermiformAppendix, Sample, ppb, colnames(VermiformAppendix)[2]:colnames(VermiformAppendix)[ncol(VermiformAppendix)], factor_key=TRUE)


#Removed Breast
All_ppb_long <- rbind(AdiposeTissue_ppb_long,AdrenalGland_ppb_long,BoneMarrow_ppb_long,Brain_ppb_long,
                       Colon_ppb_long,Duodenum_ppb_long,Esophagus_ppb_long,FallopianTubeOviduct_ppb_long,GallBladder_ppb_long,
                       Heart_ppb_long,Kidney_ppb_long,Liver_ppb_long,Lung_ppb_long,LymphNode_ppb_long,Ovary_ppb_long,
                       Pancreas_ppb_long,Placenta_ppb_long,Prostate_ppb_long,Rectum_ppb_long,SalivaryGland_ppb_long,
                       SmallIntestine_ppb_long,SmallIntestine_ppb_long,SmoothMuscle_ppb_long,Spleen_ppb_long,Stomach_ppb_long,
                       Testis_ppb_long,Thyroid_ppb_long,Tonsil_ppb_long,UmbilicalArtery_ppb_long,UrinaryBladder_ppb_long,
                       UterineEndometrium_ppb_long,VermiformAppendix_ppb_long)

All_ppb_long <- All_ppb_long[complete.cases(All_ppb_long),]
All_ppb_long$Tissues <- gsub(".*\\.", "", All_ppb_long$Sample, perl=TRUE)

All_ppb_long$Organs <- All_ppb_long$Tissues
All_ppb_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|MiddleFrontalGyrus|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|Neocortex|OccipitalCortex|PinealGland|PituitaryHypophysis|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                             "Brain", ignore.case = FALSE, All_ppb_long$Organs)
All_ppb_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                             "Heart", ignore.case = FALSE, All_ppb_long$Organs)

All_ppb_long <- merge(x=All_ppb_long, y=sample_sizes,
                       by.x=c("Organs"), by.y=c("Organs"))

All_ppb_long$Organs_samples <- paste(All_ppb_long$Organs, " (", All_ppb_long$Number_of_samples, ")", sep="")

All_ppb_long$Datasets <- gsub("\\..*", "", All_ppb_long$Sample, perl=TRUE)

ggplot(All_ppb_long, aes(x=Organs_samples, y=ppb)) + 
  #geom_violin(trim = FALSE) + 
  #geom_boxplot(width = 0.2)+
  geom_boxplot() + 
  xlab("Organs")+
  ylab("protein abundance (ppb)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("ppb")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

labs <- data.frame(vals = c(min(All_ppb_long$ppb), median(All_ppb_long$ppb), max(All_ppb_long$ppb)))
ggplot(All_ppb_long, aes(x=ppb)) + 
  geom_density() + 
  xlab("protein abundance (ppb)")+
  scale_x_log10()+
  theme_bw()+
  geom_vline(xintercept = median(All_ppb_long$ppb))+
  #geom_vline(xintercept = min(All_ppb_long$ppb))+
  #geom_vline(xintercept = max(All_ppb_long$ppb))+
  geom_label(data=labs, aes(x=1, y=0.03, label=labs[1,1]))+
  geom_label(data=labs, aes(x=labs[2,1], y=0.4, label=labs[2,1]))+
  geom_label(data=labs, aes(x=100000000, y=0.03, label=labs[3,1]))+
  ggtitle("Distribution of normalised iBAQs (as ppb) from all human datasets")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))


# Compute the median from all samples of each organ separately
Organs_median <- as.data.frame(cbind(GeneName = AdiposeTissue$GeneName, AdiposeTissue_median = apply(as.data.frame(AdiposeTissue[,-c(1)]), 1, median, na.rm = T)))
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = AdrenalGland$GeneName, AdrenalGland_median = apply(as.data.frame(AdrenalGland[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = BoneMarrow$GeneName, BoneMarrow_median = apply(as.data.frame(BoneMarrow[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Brain$GeneName, Brain_median = apply(Brain[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
#Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Breast$GeneName, Breast_median = apply(Breast[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Colon$GeneName, Colon_median = apply(Colon[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Duodenum$GeneName, Duodenum_median = apply(as.data.frame(Duodenum[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Esophagus$GeneName, Esophagus_median = apply(as.data.frame(Esophagus[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = FallopianTubeOviduct$GeneName, FallopianTubeOviduct_median = apply(as.data.frame(FallopianTubeOviduct[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = GallBladder$GeneName, GallBladder_median = apply(as.data.frame(GallBladder[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Heart$GeneName, Heart_median = apply(Heart[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Kidney$GeneName, Kidney_median = apply(as.data.frame(Kidney[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Liver$GeneName, Liver_median = apply(Liver[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Lung$GeneName, Lung_median = apply(as.data.frame(Lung[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = LymphNode$GeneName, LymphNode_median = apply(as.data.frame(LymphNode[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Ovary$GeneName, Ovary_median = apply(as.data.frame(Ovary[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Pancreas$GeneName, Pancreas_median = apply(as.data.frame(Pancreas[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Placenta$GeneName, Placenta_median = apply(as.data.frame(Placenta[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Prostate$GeneName, Prostate_median = apply(as.data.frame(Prostate[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Rectum$GeneName, Rectum_median = apply(as.data.frame(Rectum[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = SalivaryGland$GeneName, SalivaryGland_median = apply(as.data.frame(SalivaryGland[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = SmallIntestine$GeneName, SmallIntestine_median = apply(as.data.frame(SmallIntestine[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = SmoothMuscle$GeneName, SmoothMuscle_median = apply(as.data.frame(SmoothMuscle[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Spleen$GeneName, Spleen_median = apply(as.data.frame(Spleen[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Stomach$GeneName, Stomach_median = apply(as.data.frame(Stomach[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Testis$GeneName, Testis_median = apply(as.data.frame(Testis[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Thyroid$GeneName, Thyroid_median = apply(as.data.frame(Thyroid[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Tonsil$GeneName, Tonsil_median = apply(as.data.frame(Tonsil[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = UmbilicalArtery$GeneName, UmbilicalArtery_median = apply(as.data.frame(UmbilicalArtery[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = UrinaryBladder$GeneName, UrinaryBladder_median = apply(as.data.frame(UrinaryBladder[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = UterineEndometrium$GeneName, UterineEndometrium_median = apply(as.data.frame(UterineEndometrium[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = VermiformAppendix$GeneName, VermiformAppendix_median = apply(as.data.frame(VermiformAppendix[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)


Count <- data.frame(Present_in_number_of_samples = apply(Organs_median[2:ncol(Organs_median)], 1, function(x) length(which(x != "NaN")) ))
Count$Sample_percentage <- ((Count$Present_in_number_of_samples)/(ncol(Organs_median)-1)) *100
Count <- group_by(Count, Present_in_number_of_samples) %>% mutate(Total_number_of_sample_occurences = n()) %>% mutate(Gene_percent = (Total_number_of_sample_occurences / nrow(Count))*100)

Organs_median <- cbind(Organs_median, Count)

Organs_median <- merge(x=Gene_info, y=Organs_median,
                     by.x=c("GeneName"), by.y=c("GeneName"),
                     all.x=FALSE, all.y=FALSE)

plotdata <- unique(Count)

write.table(Organs_median, file = paste("Gene_distribution_in_organs-GeneNames-Median_intensities.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )


ggplot(plotdata[plotdata$Present_in_number_of_samples !=0,], aes(x=Present_in_number_of_samples, y=Gene_percent)) + 
  geom_bar(stat="identity") + 
  xlab("Present in number of organs")+
  ylab("% of identified 'canonical proteins'")+
  theme_bw()+
  ggtitle("Distribution of 'canonical proteins' across organs")

ggplot(plotdata[plotdata$Present_in_number_of_samples !=0,], aes(x=Present_in_number_of_samples, y=Total_number_of_sample_occurences)) + 
  geom_bar(stat="identity") + 
  xlab("Present in number of organs")+
  ylab("Number of canonical proteins")+
  theme_bw()+
  ggtitle("Distribution of 'canonical proteins' across organs")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))
 

###### Gene distribution per organ
#gene_counts_per_organ <- as.data.frame(sapply(Organs_median[,-c(1,2, 35:38)], function(x) sum(!is.na(x))))
gene_counts_per_organ <- as.data.frame(sapply(Organs_median[,-c(1,2, 34:37)], function(x) sum(!is.na(x))))
colnames(gene_counts_per_organ) <- "Number_of_identified_genes"
gene_counts_per_organ <- tibble::rownames_to_column(gene_counts_per_organ, "Organs")
gene_counts_per_organ$Organs <- gsub("_median", "", gene_counts_per_organ$Organs, perl=TRUE)
gene_counts_per_organ$Percentage <- (gene_counts_per_organ$Number_of_identified_genes/nrow(Organs_median))*100

gene_counts_per_organ <- merge(x=gene_counts_per_organ, y=sample_sizes,
                               by.x=c("Organs"), by.y="Organs")


gene_counts_per_organ$Organs_samples <- paste(gene_counts_per_organ$Organs, " (", gene_counts_per_organ$Number_of_samples, ")", sep="")

write.table(gene_counts_per_organ, file = paste("Gene_distribution_in_organs_plot.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

ggplot(gene_counts_per_organ, aes(x=Organs_samples, y=Percentage)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("% identified canonical proteins")+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggplot(gene_counts_per_organ, aes(x=Organs_samples, y=Number_of_identified_genes)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

##### Gene distribution per dataset
All_datasets <- Merged_input_data
#All_datasets <- aggregate(All_datasets[,-c(1,2)], list("GeneName"= All_datasets$GeneName), median, na.rm =TRUE)

dataset_sample_names <- colnames(All_datasets[-c(1)])
colnames(All_datasets) <- gsub("\\..*", "", colnames(All_datasets), perl=TRUE)

Datasets_median <- as.data.frame(cbind(GeneName = All_datasets$GeneName, PXD000547_median = apply(as.data.frame(All_datasets[,grep("PXD000547", colnames(All_datasets))]), 1, median, na.rm = T)))
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD000548_median = apply(as.data.frame(All_datasets[,grep("PXD000548", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
#Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD001325_median = apply(as.data.frame(All_datasets[,grep("PXD001325", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD001608_median = apply(as.data.frame(All_datasets[,grep("PXD001608", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD002029_median = apply(as.data.frame(All_datasets[,grep("PXD002029", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD004143_median = apply(as.data.frame(All_datasets[,grep("PXD004143", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD004332_median = apply(as.data.frame(All_datasets[,grep("PXD004332", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD005819_median = apply(as.data.frame(All_datasets[,grep("PXD005819", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD006233_median = apply(as.data.frame(All_datasets[,grep("PXD006233", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD006675_median = apply(as.data.frame(All_datasets[,grep("PXD006675", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD008934_median = apply(as.data.frame(All_datasets[,grep("PXD008934", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD010154_median = apply(as.data.frame(All_datasets[,grep("PXD010154", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD010271_median = apply(as.data.frame(All_datasets[,grep("PXD010271", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD012131_median = apply(as.data.frame(All_datasets[,grep("PXD012131", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD012755_median = apply(as.data.frame(All_datasets[,grep("PXD012755", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD015079_median = apply(as.data.frame(All_datasets[,grep("PXD015079", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD020187_median = apply(as.data.frame(All_datasets[,grep("PXD020187", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn6038852_median = apply(as.data.frame(All_datasets[,grep("syn6038852", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn21444980_median = apply(as.data.frame(All_datasets[,grep("syn21444980", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn7204174_median = apply(as.data.frame(All_datasets[,grep("syn7204174", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn3606087_median = apply(as.data.frame(All_datasets[,grep("syn3606087", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn4624471_median = apply(as.data.frame(All_datasets[,grep("syn4624471", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn7431984_median = apply(as.data.frame(All_datasets[,grep("syn7431984", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn6038797_median = apply(as.data.frame(All_datasets[,grep("syn6038797", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, syn21443008_median = apply(as.data.frame(All_datasets[,grep("syn21443008", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)


# Gene distribution per dataset
gene_counts_per_dataset <- as.data.frame(sapply(Datasets_median[,-c(1)], function(x) sum(!is.na(x))))
colnames(gene_counts_per_dataset) <- "Number_of_identified_genes"
gene_counts_per_dataset <- tibble::rownames_to_column(gene_counts_per_dataset, "Datasets")
gene_counts_per_dataset$Datasets <- gsub("_median", "", gene_counts_per_dataset$Datasets, perl=TRUE)
gene_counts_per_dataset$Percentage <- (gene_counts_per_dataset$Number_of_identified_genes/nrow(Datasets_median))*100

# Number of tissues per dataset
dataset_tissues <- data.frame(Datasets= "PXD000547", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD000547", dataset_sample_names)], perl=TRUE))) )
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD000548", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD000548", dataset_sample_names)], perl=TRUE))) ))
#dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD001325", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD001325", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD001608", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD001608", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD002029", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD002029", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD004143", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD004143", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD004332", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD004332", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD005819", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD005819", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD006233", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD006233", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD006675", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD006675", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD008934", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD008934", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD010154", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD010154", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD010271", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD010271", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD012131", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD012131", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD012755", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD012755", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD015079", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD015079", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD020187", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD020187", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn6038852", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn6038852", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn21444980", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn21444980", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn7204174", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn7204174", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn3606087", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn3606087", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn4624471", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn4624471", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn7431984", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn7431984", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn6038797", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn6038797", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "syn21443008", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("syn21443008", dataset_sample_names)], perl=TRUE))) ))


gene_counts_per_dataset <- merge(x=gene_counts_per_dataset, y=dataset_tissues,
                                 by.x=c("Datasets"), by.y=c("Datasets"), all.x=FALSE, all.y=FALSE)

gene_counts_per_dataset$Datasets_tissues <- paste(gene_counts_per_dataset$Datasets, " (", gene_counts_per_dataset$Tissue_count, ")", sep="")

write.table(gene_counts_per_dataset, file = paste("Gene_distribution_in_datasets_plot.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

ggplot(gene_counts_per_dataset, aes(x=Datasets_tissues, y=Percentage)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("% identified canonical proteins")+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across datasets")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggplot(gene_counts_per_dataset, aes(x=Datasets_tissues, y=Number_of_identified_genes)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across datasets")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

# Plot distribution of iBAQ values across datasets
All_iBAQ_long <- merge(x=All_ppb_long, y=dataset_tissues,
                       by.x=c("Datasets"), by.y=c("Datasets"))

All_iBAQ_long$Dataset_tissues <- paste(All_iBAQ_long$Datasets, " (", All_iBAQ_long$Tissue_count, ")", sep="")


ggplot(All_iBAQ_long, aes(x=Dataset_tissues, y=ppb)) + 
  #geom_violin(trim = FALSE) + 
  #geom_boxplot(width = 0.2)+
  geom_boxplot() + 
  xlab("Datasets")+
  ylab("protein abundance (ppb)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("iBAQ")+
  geom_signif(comparisons = list(c("versicolor", "virginica")), 
              map_signif_level=TRUE)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))


# Relation between sample size and gene coverage
scatterplot <- gene_counts_per_organ[,c("Percentage", "Organs_samples")]
scatterplot$Name <- gsub(" .*", "", scatterplot$Organs_samples, perl=TRUE)
scatterplot$Organs_samples <- gsub(".*\\(|\\)", "", scatterplot$Organs_samples, perl=TRUE)

  
ggplot(scatterplot, aes(x=as.numeric(Organs_samples), y=Percentage)) + 
  geom_point() + 
  geom_smooth(method='lm')+
  xlab("sample size")+
  ylab("% canonical protein coverage")+
  theme_bw()+
  scale_x_log10()+
  geom_label_repel(aes(label = Name), size = 3, max.overlaps = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Relation between sample size and 'canonical protein' coverage")+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))



Brain_long <- gather(Brain, Samples, Intensities, colnames(Brain[2]):colnames(Brain[ncol(Brain)]), factor_key=TRUE)
Brain_long$Tissues <- gsub(".*\\.", "", Brain_long$Samples, perl=TRUE)

# If gene is identified in at least one tissue sample then that gene is counted as identified
Brain_long <- Brain_long %>% group_by(GeneName, Tissues) %>% mutate(is_gene_NA_in_all_tissues = all(is.na(Intensities)))

filtered_data_long <- Brain_long[Brain_long$is_gene_NA_in_all_tissues == "FALSE",]



#Heart_long <- gather(Heart, Samples, Intensities, colnames(Heart[2]):colnames(Heart[ncol(Heart)]), factor_key=TRUE)
#Heart_long$Tissues <- gsub(".*\\.", "", Heart_long$Samples, perl=TRUE)

# If gene is identified in at least one tissue sample then that gene is counted as identified
#Heart_long <- Heart_long %>% group_by(GeneName, Tissues) %>% mutate(is_gene_NA_in_all_tissues = all(is.na(Intensities)))

#filtered_data_long <- Heart_long[Heart_long$is_gene_NA_in_all_tissues == "FALSE",]




### 
# The undetected genes are also put into bin1 
# (After discussing with Andy Jones [20/01/2021]: since undetected genes could be below detection threshold due to their low expression or abundance)
#replace missing values (NA) of those genes that were detected in at least one tissue sample
#filtered_data_long[is.na(filtered_data_long)] <- 0

filtered_data <- spread(filtered_data_long[,-c(4,5)], Samples, Intensities)

# count number of samples for each gene that have non NA values (i.e., detected)
non_missing_sample_count_percentage <- data.frame(non_missing_sample_count_percentage= apply(filtered_data[2:ncol(filtered_data)], 1, function(x) (length(which(!is.na(x))))/(ncol(filtered_data)-1)*100 ))

filtered_data <- cbind(filtered_data, non_missing_sample_count_percentage)

# Filter genes that are present in at least 33%, 50% or 75% of samples
filtered_data <- filtered_data[filtered_data$non_missing_sample_count_percentage >= 50,]


# PCA of only intensity values without batch correction
pca_input <- filtered_data[,-c(1,ncol(filtered_data))]
pca_input[is.na(pca_input)] <- 0
 


#### Normalisation 
norm_input <- filtered_data[,-c(1,ncol(filtered_data))]
norm_input[is.na(norm_input)] <- 0

# Remove tissue samples with only 1 replicate
# norm_input <- filtered_data[,-c(1,21,22,33,34,35,40,41,45,46,47, ncol(filtered_data))]

sample_names_tissues <- gsub(".*\\.", "", colnames(norm_input), perl=TRUE)
sample_names_tissues <- gsub("breast", "Breast", sample_names_tissues, perl=TRUE)
input_batch_tissues <- data.frame(sample_names_tissues)
input_batch_tissues <- input_batch_tissues %>% group_by(sample_names_tissues) %>% mutate( batchID = cur_group_id() )


sample_names_datasets <- gsub("\\..*", "", colnames(norm_input), perl=TRUE)
sample_names_datasets <- gsub("breast", "Breast", sample_names_datasets, perl=TRUE)
input_batch_datasets <- data.frame(sample_names_datasets)
input_batch_datasets <- input_batch_datasets %>% group_by(sample_names_datasets) %>% mutate( batchID = cur_group_id() )


#### Combat normalisation to remove batch effects
combat_normalised <- ComBat(as.matrix(norm_input),
                            batch = input_batch_datasets$batchID,
                            mod = NULL,
                            par.prior = TRUE,
                            prior.plots = FALSE,
                            mean.only = TRUE,
                            ref.batch = NULL,
                            BPPARAM = bpparam("SerialParam"))

# PCA of combat normalised intensity values
pca_input <- combat_normalised


####### LIMMA normalisation to remove batch effect
# https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf (page 192)


limma_normalised <-  removeBatchEffect(norm_input, 
                                       batch=input_batch_datasets$batchID, 
                                       batch2=NULL, 
                                       covariates=NULL)

# PCA of limma normalised intensity values
pca_input <- limma_normalised


#### PCA plot
# 
pca_bins <- prcomp(t(pca_input), scale = FALSE)
pca_plot_data <- data.frame(pca_bins$x[,1:2]) # Take components 1 and 2
pca_plot_data <- tibble::rownames_to_column(pca_plot_data, "Samples")

pca_plot_data$Tissues <- gsub(".*\\.", "", pca_plot_data$Samples, perl=TRUE)
pca_plot_data <- pca_plot_data %>% group_by(Tissues) %>% mutate( TissueID = cur_group_id() )

pca_plot_data$Datasets <- gsub("\\..*", "", pca_plot_data$Samples, perl=TRUE)
pca_plot_data <- pca_plot_data %>% group_by(Datasets) %>% mutate( DatasetID = cur_group_id() )

#Adding number of datasets next to each tissue sample
pca_plot_data <- pca_plot_data %>% group_by(Tissues) %>% mutate( Tissues = paste(Tissues, "(", length(unique(Datasets)), ")") )


# To add lables on the legend
# From here: https://stackoverflow.com/questions/49965758/change-geom-texts-default-a-legend-to-label-string-itself

oldK <- GeomText$draw_key
# define new key
# if you manually add colours then add vector of colours 
# instead of `scales::hue_pal()(length(var))`
GeomText$draw_key <- function (data, params, size, 
                               var=unique(pca_plot_data$DatasetID), 
                               cols=scales::hue_pal()(length(var))) {
  
  # sort as ggplot sorts these alphanumerically / or levels of factor
  txt <- if(is.factor(var)) levels(var) else sort(var)
  txt <- txt[match(data$colour, cols)]
  
  textGrob(txt, 0.5, 0.5,  
           just="center", 
           gp = gpar(col = alpha(data$colour, data$alpha), 
                     fontfamily = data$family, 
                     fontface = data$fontface, 
                     fontsize = data$size * .pt))
}

ggplot(pca_plot_data, aes(x=PC1, y = PC2, colour = Datasets))+
  #geom_point(alpha=0.6)+
  geom_text(aes(label=DatasetID))+
  labs(x="PC1", y="PC2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #labs(color="Samples") +
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.3,"line"))+
  guides(col = guide_legend(ncol = 3))+
  ggtitle("Heart-Samples_binned-by-regions batch-per-datasets\n[filter: genes detected in at least 50% of samples]")

# IMPORTANT reset key
GeomText$draw_key <- oldK


#foo <- Merged_input_data[,c(1:6)] %>%
#       group_by(GeneID, GeneName) %>% 
#       mutate(PXD000547 = all(is.na(c(PXD000547.Sample1.CorpusCallosum, PXD000547.Sample2.CorpusCallosum))), 
#              PXD000548 = all(is.na(c(PXD000548.Sample1.AnteriorTemporalLobe, PXD000548.Sample2.AnteriorTemporalLobe))))
