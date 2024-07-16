library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(ggpubr)

##### CRC Secretome data #########
MS_spectra_data_secretome <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Secretome/All_DDA_datasets_Spectra_counts_secretome.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
sample_classification_secretome <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Secretome/Sample_classification.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
sample_classification_secretome$Experiment <- gsub("ppb.iBAQ.","",sample_classification_secretome$Sample, perl=TRUE)
sample_classification_secretome$Experiment <- gsub("\\."," ",sample_classification_secretome$Experiment, perl=TRUE)
sample_classification_secretome$Experiment <- gsub("e EV","e-EV",sample_classification_secretome$Experiment, perl=TRUE)
sample_classification_secretome$Experiment <- gsub("CRC LM","CRC-LM",sample_classification_secretome$Experiment, perl=TRUE)

merged_secretome <- merge(x=MS_spectra_data_secretome, y=sample_classification_secretome,
                      by.x=c("Experiment"), by.y=c("Experiment"),
                      all.x=TRUE, all.y=TRUE)

merged_secretome <- merged_secretome[-1,]

secretome_MS_spectra_tissue <- merged_secretome[,c("MS","MS.MS","MS.MS.identified","Source")]
Spectra_sum_secretome_tissue <- aggregate(secretome_MS_spectra_tissue[,c(1:3)], by=list(Tissue=secretome_MS_spectra_tissue$Source), FUN=mean)
colnames(Spectra_sum_secretome_tissue) <- c("Tissues", "number of MS spectra recorded", "number of MS/MS spectra recorded", "total number of identified\ntandem MS spectra")

secretome_MS_spectra_dataset <- merged_secretome[,c("MS","MS.MS","MS.MS.identified","Dataset")]
Spectra_sum_secretome_dataset <- aggregate(secretome_MS_spectra_dataset[,c(1:3)], by=list(Dataset=secretome_MS_spectra_dataset$Dataset), FUN=mean)
colnames(Spectra_sum_secretome_dataset) <- c("Dataset", "number of MS spectra recorded", "number of MS/MS spectra recorded", "total number of identified\ntandem MS spectra")

Canonical_proteins_secretome_tissue <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Secretome/Gene_distribution_in_tissues_plot-CRC_Secretome.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
Canonical_proteins_secretome_tissue$Tissues <- gsub("Blood_derived", "Blood-derived",Canonical_proteins_secretome_tissue$Tissues)
Canonical_proteins_secretome_tissue$Tissues <- gsub("Cellculture_derived", "Cell culture-derived",Canonical_proteins_secretome_tissue$Tissues)
Canonical_proteins_secretome_tissue$Tissues <- gsub("ECvesicles_and_exosomes", "EC vesicles and exosomes",Canonical_proteins_secretome_tissue$Tissues)
Canonical_proteins_secretome_tissue$Tissues <- gsub("Interstitial_fluid", "Interstitial fluid",Canonical_proteins_secretome_tissue$Tissues)

Spectra_sum_secretome_tissue <- merge(x=Spectra_sum_secretome_tissue, y=Canonical_proteins_secretome_tissue,
                                  by.x=c("Tissues"), by.y=c("Tissues"))


Canonical_proteins_secretome_dataset <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Secretome/Gene_distribution_in_datasets_plot-CRC_Secretome.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)

Spectra_sum_secretome_dataset <- merge(x=Spectra_sum_secretome_dataset, y=Canonical_proteins_secretome_dataset,
                                   by.x=c("Dataset"), by.y=c("Datasets"))

Spectra_sum_secretome_tissue_long <- gather(Spectra_sum_secretome_tissue, Index, Value, colnames(Spectra_sum_secretome_tissue)[2:4], factor_key=TRUE)
Spectra_sum_secretome_tissue_long$Type <- rep("Secretomes", nrow(Spectra_sum_secretome_tissue_long))
Spectra_sum_secretome_tissue_long_subset <- Spectra_sum_secretome_tissue_long[,c(2,1,6:8)]
colnames(Spectra_sum_secretome_tissue_long_subset)[2] <- "Samples"

Spectra_sum_secretome_dataset_long <- gather(Spectra_sum_secretome_dataset, Index, Value, colnames(Spectra_sum_secretome_dataset)[2:4], factor_key=TRUE)
Spectra_sum_secretome_dataset_long$Type <- rep("Datasets", nrow(Spectra_sum_secretome_dataset_long))
Spectra_sum_secretome_dataset_long_subset <- Spectra_sum_secretome_dataset_long[,c(2,1,6:8)]
colnames(Spectra_sum_secretome_dataset_long_subset)[2] <- "Samples"

plotdata_secretome <- rbind(Spectra_sum_secretome_tissue_long_subset, Spectra_sum_secretome_dataset_long_subset)

ggplot(plotdata_secretome, aes(x=Value, y=Number_of_identified_genes)) + 
  geom_point() + 
  xlab("Total spectra counts")+
  ylab("Identified canonical proteins")+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_label_repel(aes(label = Samples), size = 3, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Identified proteins vs. Quantity of data")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  #  facet_grid(Type~Index)+
  facet_grid(Type~Index, scales = "free_y")

##### CRC Tumor data #########
MS_spectra_data_tumor <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Tumour/All_DDA_datasets_Spectra_counts_tumor.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
sample_classification_tumor <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Tumour/Sample_classification.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
sample_classification_tumor$Experiment <- gsub("ppb.iBAQ.","",sample_classification_tumor$Name, perl=TRUE)
sample_classification_tumor$Experiment <- gsub("Non.tumor.","Non-tumor ",sample_classification_tumor$Experiment, perl=TRUE)
sample_classification_tumor$Experiment <- gsub("Tumor.","Tumor ",sample_classification_tumor$Experiment, perl=TRUE)
sample_classification_tumor$Experiment <- gsub("\\.","-",sample_classification_tumor$Experiment, perl=TRUE)

merged_tumor <- merge(x=MS_spectra_data_tumor, y=sample_classification_tumor,
                   by.x=c("Experiment"), by.y=c("Experiment"),
                   all.x=TRUE, all.y=TRUE)

merged_tumor <- merged_tumor[-1,]
                                     
tumor_MS_spectra_tissue <- merged_tumor[,c("MS","MS.MS","MS.MS.identified","Tissue")]
Spectra_sum_Tumor_tissue <- aggregate(tumor_MS_spectra_tissue[,c(1:3)], by=list(Tissue=tumor_MS_spectra_tissue$Tissue), FUN=mean)
colnames(Spectra_sum_Tumor_tissue) <- c("Tissues", "number of MS spectra recorded", "number of MS/MS spectra recorded", "total number of identified\ntandem MS spectra")
Spectra_sum_Tumor_tissue$Tissues <- gsub("Adenoma", "Colorectal_adenoma", Spectra_sum_Tumor_tissue$Tissues)
Spectra_sum_Tumor_tissue$Tissues <- gsub("Tumor", "Colorectal_tumor", Spectra_sum_Tumor_tissue$Tissues)
Spectra_sum_Tumor_tissue$Tissues <- gsub("Mucosa", "Colorectal_mucosa", Spectra_sum_Tumor_tissue$Tissues)


tumor_MS_spectra_dataset <- merged_tumor[,c("MS","MS.MS","MS.MS.identified","Dataset")]
Spectra_sum_Tumor_dataset <- aggregate(tumor_MS_spectra_dataset[,c(1:3)], by=list(Dataset=tumor_MS_spectra_dataset$Dataset), FUN=mean)
colnames(Spectra_sum_Tumor_dataset) <- c("Dataset", "number of MS spectra recorded", "number of MS/MS spectra recorded", "total number of identified\ntandem MS spectra")

Canonical_proteins_tumor_tissue <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Tumour/Gene_distribution_in_tissues_plot-CRC_Tumor.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
Spectra_sum_Tumor_tissue <- merge(x=Spectra_sum_Tumor_tissue, y=Canonical_proteins_tumor_tissue,
                                  by.x=c("Tissues"), by.y=c("Tissues"))

Canonical_proteins_tumor_dataset <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Tumour/Gene_distribution_in_datasets_plot-CRC_Tumor.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
Spectra_sum_Tumor_dataset <- merge(x=Spectra_sum_Tumor_dataset, y=Canonical_proteins_tumor_dataset,
                                  by.x=c("Dataset"), by.y=c("Datasets"))


Spectra_sum_Tumor_tissue_long <- gather(Spectra_sum_Tumor_tissue, Index, Value, colnames(Spectra_sum_Tumor_tissue)[2:4], factor_key=TRUE)
Spectra_sum_Tumor_tissue_long$Type <- rep("Tissues", nrow(Spectra_sum_Tumor_tissue_long))
Spectra_sum_Tumor_tissue_long_subset <- Spectra_sum_Tumor_tissue_long[,c(2,1,6:8)]
colnames(Spectra_sum_Tumor_tissue_long_subset)[2] <- "Samples"

Spectra_sum_Tumor_dataset_long <- gather(Spectra_sum_Tumor_dataset, Index, Value, colnames(Spectra_sum_Tumor_dataset)[2:4], factor_key=TRUE)
Spectra_sum_Tumor_dataset_long$Type <- rep("Datasets", nrow(Spectra_sum_Tumor_dataset_long))
Spectra_sum_Tumor_dataset_long_subset <- Spectra_sum_Tumor_dataset_long[,c(2,1,6:8)]
colnames(Spectra_sum_Tumor_dataset_long_subset)[2] <- "Samples"

plotdata_tumor <- rbind(Spectra_sum_Tumor_tissue_long_subset, Spectra_sum_Tumor_dataset_long_subset)
plotdata_tumor$Samples <- gsub("_"," ", plotdata_tumor$Samples)

ggplot(plotdata_tumor, aes(x=Value, y=Number_of_identified_genes)) + 
  geom_point() + 
  xlab("Total spectra counts")+
  ylab("Identified canonical proteins")+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_label_repel(aes(label = Samples), size = 3, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Identified proteins vs. Quantity of data")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
#  facet_grid(Type~Index)+
facet_grid(Type~Index, scales = "free_y")
