library(dplyr)
library(tidyr)
library(qdap)

MountSinai <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/proteinGroups.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

annot <- read.csv(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/ABRF2015_MaxQuant_annotation.csv", quote = "\"", header = TRUE, sep = ",", stringsAsFactors = FALSE, comment.char = "#")

annot_exclude <- annot[annot$Condition == "Exclude", c("individual")]

filtered_dataset <- MountSinai %>% select(-contains(annot_exclude))

annot$Condition <- gsub("AsymAD", "AsymptomaticAlzheimer'sDisease", annot$Condition)
annot$Condition <- gsub("AD", "Alzheimer'sDisease", annot$Condition)
annot$Samplenames <- paste("FrontalPole", annot$Condition, annot$individual, sep=".")

excluded <- annot[annot$Condition != "Exclude", ]

colnames(filtered_dataset) <- mgsub(excluded$individual, excluded$Samplenames, colnames(filtered_dataset))

write.table(filtered_dataset, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/proteinGroups.txt", sep = "\t", row.names = FALSE, quote = FALSE )

