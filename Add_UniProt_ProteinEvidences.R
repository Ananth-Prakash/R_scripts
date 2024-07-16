library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(scales)

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/")

Human_file_list <- data.frame(Dataset = c(
  # "PXD010154_Ananth", 
  # "PXD005819_33threads_yoda", 
  # "PXD004143",
  # "PXD006233",
  # "PXD012755", 
  # "PXD001608_30threads_yoda", 
  # "PXD002029", 
  # "PXD000547", 
  # "PXD000548", 
  # "PXD010271", 
  # "PXD004332",
  # "PXD006675",
  # "PXD008934", 
  # "Synapse-AD/ACT_DorsoLateralPreFrontalCortex", 
  # "Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex",
  # "Synapse-AD/Banner_DorsoLateralPreFrontalCortex",
  # "Synapse-AD/BLSA_DorsoLateralPreFrontalCortex",
  # "Synapse-AD/BLSA_Precuneus", 
  # "Synapse-AD/Mayo_TemporalCortex", 
  # "Synapse-AD/MountSinai_FrontalPole",
  # "Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases",
  # "PXD012131",
  # "PXD020187", 
  # "PXD015079"
))

Gene_ProteinID_alldatasets_bins <- read.table(file="/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/JPR/SupplementaryTable_3.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
Gene_ProteinID_alldatasets_bins <- Gene_ProteinID_alldatasets_bins[Gene_ProteinID_alldatasets_bins$UniProt.ID != "",]

## Using UniProt Protein Existence Evidences (release 2022_02)
#ProteinEvidences <- read.table(file="/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Uniprot-Human_ProteinEvidences.txt", header = TRUE, quote = "\"", sep = "\t", stringsAsFactors = FALSE)
#ProteinEvidences$Protein.existence <- gsub("Evidence at protein level","1", ProteinEvidences$Protein.existence)
#ProteinEvidences$Protein.existence <- gsub("Evidence at transcript level","2", ProteinEvidences$Protein.existence)
#ProteinEvidences$Protein.existence <- gsub("Inferred from homology","3", ProteinEvidences$Protein.existence)
#ProteinEvidences$Protein.existence <- gsub("Predicted","4", ProteinEvidences$Protein.existence)
#ProteinEvidences$Protein.existence <- gsub("Uncertain","5", ProteinEvidences$Protein.existence)


## Using NextProt Protein Existence Evidences (downloaded Oct 2022)
ProteinEvidences <- read.table(file="/Users/ananth/Documents/MaxQuant_Bechmarking/Human/NextProt-Human_ProteinEvidences.txt", header = TRUE, quote = "\"", sep = "\t", stringsAsFactors = FALSE)

merge_IDs <- function(dataset){
  mergeddata <- merge(x=dataset, y=Gene_ProteinID_alldatasets_bins[,c("Gene.Name","Gene.ID","UniProt.ID")],
                      by.x=c("Gene.Symbol","Gene.ID"), by.y=c("Gene.Name","Gene.ID"),
                      all.x=FALSE, all.y=FALSE)
}

merge_evidences <- function(df){
  merged_PE <- merge(x=df, y=ProteinEvidences[,c("Entry","Protein.existence")],
                     by.x=c("UniProt.ID"), by.y=c("Entry"),
                     all.x=TRUE, all.y=FALSE) 
}


dataset1 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010154_Ananth/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset1 <- dataset1[,c("Gene.ID","Gene.Symbol")]
dataset1 <- merge_IDs(dataset1)
dataset1 <- merge_evidences(dataset1)
colnames(dataset1)[4] <- c("ProteinEvidence")
dataset1$PX1 <- rep("PXD010154", nrow(dataset1))

dataset2 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD005819_33threads_yoda/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset2 <- dataset2[,c("Gene.ID","Gene.Name")]
colnames(dataset2)[2] <- c("Gene.Symbol")
dataset2 <- merge_IDs(dataset2)
dataset2 <- merge_evidences(dataset2)
colnames(dataset2)[4] <- c("ProteinEvidence")
dataset2$PX2 <- rep("PXD005819", nrow(dataset2))

dataset3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004143/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset3 <- dataset3[,c("Gene.ID","Gene.Name")]
colnames(dataset3)[2] <- c("Gene.Symbol")
dataset3 <- merge_IDs(dataset3)
dataset3 <- merge_evidences(dataset3)
colnames(dataset3)[4] <- c("ProteinEvidence")
dataset3$PX3 <- rep("PXD004143", nrow(dataset3))

dataset4 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006233/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset4 <- dataset4[,c("Gene.ID","Gene.Name")]
colnames(dataset4)[2] <- c("Gene.Symbol")
dataset4 <- merge_IDs(dataset4)
dataset4 <- merge_evidences(dataset4)
colnames(dataset4)[4] <- c("ProteinEvidence")
dataset4$PX4 <- rep("PXD006233", nrow(dataset4))

dataset5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012755/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset5 <- dataset5[,c("Gene.ID","Gene.Name")]
colnames(dataset5)[2] <- c("Gene.Symbol")
dataset5 <- merge_IDs(dataset5)
dataset5 <- merge_evidences(dataset5)
colnames(dataset5)[4] <- c("ProteinEvidence")
dataset5$PX5 <- rep("PXD012755", nrow(dataset5))

dataset6 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001608_30threads_yoda/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset6 <- dataset6[,c("Gene.ID","Gene.Name")]
colnames(dataset6)[2] <- c("Gene.Symbol")
dataset6 <- merge_IDs(dataset6)
dataset6 <- merge_evidences(dataset6)
colnames(dataset6)[4] <- c("ProteinEvidence")
dataset6$PX6 <- rep("PXD001608", nrow(dataset6))

dataset7 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD002029/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset7 <- dataset7[,c("Gene.ID","Gene.Name")]
colnames(dataset7)[2] <- c("Gene.Symbol")
dataset7 <- merge_IDs(dataset7)
dataset7 <- merge_evidences(dataset7)
colnames(dataset7)[4] <- c("ProteinEvidence")
dataset7$PX7 <- rep("PXD002029", nrow(dataset7))

dataset8 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000547/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset8 <- dataset8[,c("Gene.ID","Gene.Name")]
colnames(dataset8)[2] <- c("Gene.Symbol")
dataset8 <- merge_IDs(dataset8)
dataset8 <- merge_evidences(dataset8)
colnames(dataset8)[4] <- c("ProteinEvidence")
dataset8$PX8 <- rep("PXD000547", nrow(dataset8))

dataset9 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000548/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset9 <- dataset9[,c("Gene.ID","Gene.Name")]
colnames(dataset9)[2] <- c("Gene.Symbol")
dataset9 <- merge_IDs(dataset9)
dataset9 <- merge_evidences(dataset9)
colnames(dataset9)[4] <- c("ProteinEvidence")
dataset9$PX9 <- rep("PXD000548", nrow(dataset9))

dataset10 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010271/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset10 <- dataset10[,c("Gene.ID","Gene.Name")]
colnames(dataset10)[2] <- c("Gene.Symbol")
dataset10 <- merge_IDs(dataset10)
dataset10 <- merge_evidences(dataset10)
colnames(dataset10)[4] <- c("ProteinEvidence")
dataset10$PX10 <- rep("PXD010271", nrow(dataset10))

dataset11 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004332/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset11 <- dataset11[,c("Gene.ID","Gene.Name")]
colnames(dataset11)[2] <- c("Gene.Symbol")
dataset11 <- merge_IDs(dataset11)
dataset11 <- merge_evidences(dataset11)
colnames(dataset11)[4] <- c("ProteinEvidence")
dataset11$PX11 <- rep("PXD004332", nrow(dataset11))

dataset12 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006675/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset12 <- dataset12[,c("Gene.ID","Gene.Name")]
colnames(dataset12)[2] <- c("Gene.Symbol")
dataset12 <- merge_IDs(dataset12)
dataset12 <- merge_evidences(dataset12)
colnames(dataset12)[4] <- c("ProteinEvidence")
dataset12$PX12 <- rep("PXD006675", nrow(dataset12))

dataset13 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008934/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset13 <- dataset13[,c("Gene.ID","Gene.Name")]
colnames(dataset13)[2] <- c("Gene.Symbol")
dataset13 <- merge_IDs(dataset13)
dataset13 <- merge_evidences(dataset13)
colnames(dataset13)[4] <- c("ProteinEvidence")
dataset13$PX13 <- rep("PXD008934", nrow(dataset13))

dataset14 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/ACT_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset14 <- dataset14[,c("Gene.ID","Gene.Symbol")]
dataset14 <- merge_IDs(dataset14)
dataset14 <- merge_evidences(dataset14)
colnames(dataset14)[4] <- c("ProteinEvidence")
dataset14$PX14 <- rep("syn6038852", nrow(dataset14))

dataset15 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset15 <- dataset15[,c("Gene.ID","Gene.Symbol")]
dataset15 <- merge_IDs(dataset15)
dataset15 <- merge_evidences(dataset15)
colnames(dataset15)[4] <- c("ProteinEvidence")
dataset15$PX15 <- rep("syn21444980", nrow(dataset15))

dataset16 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Banner_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset16 <- dataset16[,c("Gene.ID","Gene.Symbol")]
dataset16 <- merge_IDs(dataset16)
dataset16 <- merge_evidences(dataset16)
colnames(dataset16)[4] <- c("ProteinEvidence")
dataset16$PX16 <- rep("syn7204174", nrow(dataset16))

dataset17 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset17 <- dataset17[,c("Gene.ID","Gene.Symbol")]
dataset17 <- merge_IDs(dataset17)
dataset17 <- merge_evidences(dataset17)
colnames(dataset17)[4] <- c("ProteinEvidence")
dataset17$PX17 <- rep("syn3606087", nrow(dataset17))

dataset18 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_Precuneus/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset18 <- dataset18[,c("Gene.ID","Gene.Symbol")]
dataset18 <- merge_IDs(dataset18)
dataset18 <- merge_evidences(dataset18)
colnames(dataset18)[4] <- c("ProteinEvidence")
dataset18$PX18 <- rep("syn4624471", nrow(dataset18))

dataset19 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Mayo_TemporalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset19 <- dataset19[,c("Gene.ID","Gene.Symbol")]
dataset19 <- merge_IDs(dataset19)
dataset19 <- merge_evidences(dataset19)
colnames(dataset19)[4] <- c("ProteinEvidence")
dataset19$PX19 <- rep("syn7431984", nrow(dataset19))

dataset20 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset20 <- dataset20[,c("Gene.ID","Gene.Name")]
colnames(dataset20)[2] <- c("Gene.Symbol")
dataset20 <- merge_IDs(dataset20)
dataset20 <- merge_evidences(dataset20)
colnames(dataset20)[4] <- c("ProteinEvidence")
dataset20$PX20 <- rep("syn6038797", nrow(dataset20))

dataset21 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset21 <- dataset21[,c("Gene.ID","Gene.Symbol")]
dataset21 <- merge_IDs(dataset21)
dataset21 <- merge_evidences(dataset21)
colnames(dataset21)[4] <- c("ProteinEvidence")
dataset21$PX21 <- rep("syn21443008", nrow(dataset21))

dataset22 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012131/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset22 <- dataset22[,c("Gene.ID","Gene.Symbol")]
dataset22 <- merge_IDs(dataset22)
dataset22 <- merge_evidences(dataset22)
colnames(dataset22)[4] <- c("ProteinEvidence")
dataset22$PX22 <- rep("PXD012131", nrow(dataset22))

dataset23 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD020187/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset23 <- dataset23[,c("Gene.ID","Gene.Symbol")]
dataset23 <- merge_IDs(dataset23)
dataset23 <- merge_evidences(dataset23)
colnames(dataset23)[4] <- c("ProteinEvidence")
dataset23$PX23 <- rep("PXD020187", nrow(dataset23))

dataset24 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD015079/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset24 <- dataset24[,c("Gene.ID","Gene.Symbol")]
dataset24 <- merge_IDs(dataset24)
dataset24 <- merge_evidences(dataset24)
colnames(dataset24)[4] <- c("ProteinEvidence")
dataset24$PX24 <- rep("PXD015079", nrow(dataset24))



merge_datasets <- function(df1, df2){
  data <- merge(x=df1, y=df2,
                by.x=c("UniProt.ID", "Gene.Symbol", "Gene.ID", "ProteinEvidence"),
                by.y=c("UniProt.ID", "Gene.Symbol", "Gene.ID", "ProteinEvidence"),
                all.x=TRUE, all.y=TRUE)
}

All_data <- merge_datasets(dataset1, dataset2)
All_data <- merge_datasets(All_data, dataset3)
All_data <- merge_datasets(All_data, dataset4)
All_data <- merge_datasets(All_data, dataset5)
All_data <- merge_datasets(All_data, dataset6)
All_data <- merge_datasets(All_data, dataset7)
All_data <- merge_datasets(All_data, dataset8)
All_data <- merge_datasets(All_data, dataset9)
All_data <- merge_datasets(All_data, dataset10)
All_data <- merge_datasets(All_data, dataset11)
All_data <- merge_datasets(All_data, dataset12)
All_data <- merge_datasets(All_data, dataset13)
All_data <- merge_datasets(All_data, dataset14)
All_data <- merge_datasets(All_data, dataset15)
All_data <- merge_datasets(All_data, dataset16)
All_data <- merge_datasets(All_data, dataset17)
All_data <- merge_datasets(All_data, dataset18)
All_data <- merge_datasets(All_data, dataset19)
All_data <- merge_datasets(All_data, dataset20)
All_data <- merge_datasets(All_data, dataset21)
All_data <- merge_datasets(All_data, dataset22)
All_data <- merge_datasets(All_data, dataset23)
All_data <- merge_datasets(All_data, dataset24)

All_data$Datasets <- apply(All_data[ ,c("PX1","PX2","PX3", "PX4","PX5","PX6","PX7",
                                        "PX8","PX9","PX10","PX11","PX12","PX13","PX14",
                                        "PX15","PX16","PX17","PX18","PX19","PX20","PX21",
                                        "PX22","PX23","PX24") ] , 1 , paste , collapse = "," )
All_data$Present_in_number_of_datasets <- apply(All_data[5:28], 1, function(x) length(which(!is.na(x))) )

All_data$Datasets <- gsub("NA,","", All_data$Datasets)
All_data$Datasets <- gsub(",NA$","",All_data$Datasets, perl=TRUE)

All_data <- All_data[,-c(5:28)]

#Some of the UniProt entries have been made obsolete in the latest release 2022_02, 
#therefore manually adding Protein Evidence values
#All_data$UniProt_ProteinEvidence[All_data$UniProt.ID == "A0A499FJT0"] <- 1
#All_data$UniProt_ProteinEvidence[All_data$UniProt.ID == "H0YGT0"] <- 4
#All_data$UniProt_ProteinEvidence[All_data$UniProt.ID == "Q59FZ7"] <- 2
#All_data$UniProt_ProteinEvidence[All_data$UniProt.ID == "Q6ZW33"] <- 1

write.table(All_data, file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/NextProt_ProteinEvidencesDistribution_per_dataset.txt", sep = "\t", row.names = FALSE, quote = FALSE )


foo <- All_data[complete.cases(All_data),]
foo <- foo[foo$UniProt.ID != "B1AH88",]

ggplot(foo, aes(fill=as.factor(ProteinEvidence), x=Present_in_number_of_datasets)) + 
  geom_bar(colour="grey")+ 
  xlab("Number of datasets")+
  ylab("Number of identified canonical proteins")+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(fill="NextProt Protein Existence")+
  scale_fill_manual(values = c("#F8766D", "#A3A500", "#00B0F6"))+
  theme(legend.position = c(0.8, 0.8))+
  ggtitle("NextProt Protein Existence Evidence distribution in datasets")


ggplot(All_data[All_data$ProteinEvidence > 1, ], aes(fill=as.factor(ProteinEvidence), x=Present_in_number_of_datasets)) + 
  geom_bar(colour="grey")+ 
  xlab("Number of datasets")+
  ylab("Number of identified canonical proteins")+
  #scale_y_log10()+
  theme_bw()+
  labs(fill="NextProt Protein Existence")+
  scale_fill_manual(values = c("#A3A500", "#00B0F6"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = c(0.8, 0.8))+
  ggtitle("NextProt Protein Existence Evidence distribution in datasets\n(PE > 1)")



