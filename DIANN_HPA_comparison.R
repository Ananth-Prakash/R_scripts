library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(gplots)

# Comparison with Human Protein Atlas (HPA)
# Data download: https://www.proteinatlas.org/about/download
# 1	Normal tissue data. Date downloaded: 6/Jan/2022

# Read HPA normal tissue protein expression based on 
# "Expression profiles for proteins in human tissues based on immunohistochemistry using tissue micro arrays"

setwd("/Users/ananth/Documents/DIANN/")


tmp  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/normal_tissue.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# details of Reliability scores are here: 
# https://www.proteinatlas.org/about/antibody+validation#enhanced_antibody_validation___ihc

# remove entries with Reliabilty == Uncertain

tmp_filtered <- tmp[tmp$Reliability != "Uncertain",]

# Replace 'Level' annotations with scores that are similar to bin values, useful for aggregating over cell types
# "Not detected" <- NA
# "Low" <- 1
# "Medium" <- 3 
# "High" <- 5
# "Ascending" <- 1          (171 entries)
# "Descending" <- 1         (73 entries)
# "Not representative" <- 1 (23 entries)

tmp_filtered$Level[tmp_filtered$Level == "Not detected"] <- 0
tmp_filtered$Level[tmp_filtered$Level == "Low"] <- 1
tmp_filtered$Level[tmp_filtered$Level == "Medium"] <- 2
tmp_filtered$Level[tmp_filtered$Level == "High"] <- 3
tmp_filtered$Level[tmp_filtered$Level == "Ascending"] <- 1
tmp_filtered$Level[tmp_filtered$Level == "Descending"] <- 1
tmp_filtered$Level[tmp_filtered$Level == "Not representative"] <- 1

tmp_filtered$Level <- as.numeric(tmp_filtered$Level)
tmp_filtered$Level[tmp_filtered$Level == 0] <- NA

TissueTypes <- data.frame("Tissue" = unique(tmp_filtered$Tissue))

tmp_filtered$Organ <- tmp_filtered$Tissue
tmp_filtered$Organ <- gsub("caudate|cerebellum|cerebral cortex|hippocampus|hypothalamus|pituitary gland|dorsal raphe|choroid plexus|substantia nigra", "Brain", tmp_filtered$Organ, perl=TRUE)
tmp_filtered$Organ <- gsub("endometrium 1|endometrium 2", "Endometrium", tmp_filtered$Organ, perl=TRUE)
tmp_filtered$Organ <- gsub("skin 1|skin 2|skin", "Skin", tmp_filtered$Organ, perl=TRUE)
tmp_filtered$Organ <- gsub("soft tissue 1|soft tissue 2", "Soft tissue", tmp_filtered$Organ, perl=TRUE)
tmp_filtered$Organ <- gsub("stomach 1|stomach 2", "Stomach", tmp_filtered$Organ, perl=TRUE)
tmp_filtered$Organ <- gsub("retina|eye", "Eye", tmp_filtered$Organ, perl=TRUE)

OrganTypes <- data.frame("Organs" = unique(tmp_filtered$Organ))

tmp_subdata <- tmp_filtered[,c("Gene","Gene.name","Level","Organ")]

tmp_subdata_aggregate <- aggregate(tmp_subdata$Level, by = list(tmp_subdata$Gene, tmp_subdata$Gene.name, tmp_subdata$Organ), FUN = median, na.rm =TRUE)

colnames(tmp_subdata_aggregate) <- c("Gene","Gene.name","Organ","Level")

HPA_median_bins <- spread(tmp_subdata_aggregate, Organ, Level)

colnames(HPA_median_bins) <- gsub("$", "_HPA", colnames(HPA_median_bins), perl=TRUE)
colnames(HPA_median_bins) <- gsub(" ", "_", colnames(HPA_median_bins), perl=TRUE)

DIA_3bins <- read.table( "/Users/ananth/Documents/DIANN/report_AllDatasets_ControlSamplesOnly_3binsOnly_HPAcomparison.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

Gene_info <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-GeneNames-Median_intensities.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
Gene_info <- Gene_info[,c(1,2)]

DIA_3bins_gene <- merge(x=Gene_info, y=DIA_3bins,
                           by.x=c("GeneName"), by.y=c("Genes"))

colnames(DIA_3bins_gene) <- gsub("$", "_DIA", colnames(DIA_3bins_gene), perl=TRUE)

DIA_Brain <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Brain_DIA),c("GeneID_DIA","GeneName_DIA","Brain_DIA"), drop=FALSE]
HPA_Brain <- HPA_median_bins[!is.na(HPA_median_bins$Brain_HPA),c("Gene_HPA","Gene.name_HPA","Brain_HPA"), drop=FALSE]
colnames(DIA_Brain) <- c("GeneID","GeneName","Level")
colnames(HPA_Brain) <- c("GeneID","GeneName","Level")
Total_Brain <- unique(rbind(DIA_Brain[,c("GeneID"), drop=FALSE], HPA_Brain[,c("GeneID"), drop=FALSE]))
Common_Brain        <- subset(DIA_Brain, (GeneID %in% HPA_Brain$GeneID))
Only_in_DIA_Brain   <- subset(DIA_Brain, !(GeneID %in% HPA_Brain$GeneID))
Only_in_HPA_Brain   <- subset(HPA_Brain, !(GeneID %in% DIA_Brain$GeneID))

DIA_Colon <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Colon_DIA),c("GeneID_DIA","GeneName_DIA","Colon_DIA"), drop=FALSE]
HPA_Colon <- HPA_median_bins[!is.na(HPA_median_bins$colon_HPA),c("Gene_HPA","Gene.name_HPA","colon_HPA"), drop=FALSE]
colnames(DIA_Colon) <- c("GeneID","GeneName","Level")
colnames(HPA_Colon) <- c("GeneID","GeneName","Level")
Total_Colon <- unique(rbind(DIA_Colon[,c("GeneID"), drop=FALSE], HPA_Colon[,c("GeneID"), drop=FALSE]))
Common_Colon        <- subset(DIA_Colon, (GeneID %in% HPA_Colon$GeneID))
Only_in_DIA_Colon   <- subset(DIA_Colon, !(GeneID %in% HPA_Colon$GeneID))
Only_in_HPA_Colon   <- subset(HPA_Colon, !(GeneID %in% DIA_Colon$GeneID))

DIA_Duodenum <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Duodenum_DIA),c("GeneID_DIA","GeneName_DIA","Duodenum_DIA"), drop=FALSE]
HPA_Duodenum <- HPA_median_bins[!is.na(HPA_median_bins$duodenum_HPA),c("Gene_HPA","Gene.name_HPA","duodenum_HPA"), drop=FALSE]
colnames(DIA_Duodenum) <- c("GeneID","GeneName","Level")
colnames(HPA_Duodenum) <- c("GeneID","GeneName","Level")
Total_Duodenum <- unique(rbind(DIA_Duodenum[,c("GeneID"), drop=FALSE], HPA_Duodenum[,c("GeneID"), drop=FALSE]))
Common_Duodenum        <- subset(DIA_Duodenum, (GeneID %in% HPA_Duodenum$GeneID))
Only_in_DIA_Duodenum   <- subset(DIA_Duodenum, !(GeneID %in% HPA_Duodenum$GeneID))
Only_in_HPA_Duodenum   <- subset(HPA_Duodenum, !(GeneID %in% DIA_Duodenum$GeneID))

DIA_Esophagus <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Esophagus_DIA),c("GeneID_DIA","GeneName_DIA","Esophagus_DIA"), drop=FALSE]
HPA_Esophagus <- HPA_median_bins[!is.na(HPA_median_bins$esophagus_HPA),c("Gene_HPA","Gene.name_HPA","esophagus_HPA"), drop=FALSE]
colnames(DIA_Esophagus) <- c("GeneID","GeneName","Level")
colnames(HPA_Esophagus) <- c("GeneID","GeneName","Level")
Total_Esophagus <- unique(rbind(DIA_Esophagus[,c("GeneID"), drop=FALSE], HPA_Esophagus[,c("GeneID"), drop=FALSE]))
Common_Esophagus        <- subset(DIA_Esophagus, (GeneID %in% HPA_Esophagus$GeneID))
Only_in_DIA_Esophagus   <- subset(DIA_Esophagus, !(GeneID %in% HPA_Esophagus$GeneID))
Only_in_HPA_Esophagus   <- subset(HPA_Esophagus, !(GeneID %in% DIA_Esophagus$GeneID))

DIA_Heart <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Heart_DIA),c("GeneID_DIA","GeneName_DIA","Heart_DIA"), drop=FALSE]
HPA_Heart <- HPA_median_bins[!is.na(HPA_median_bins$heart_muscle_HPA),c("Gene_HPA","Gene.name_HPA","heart_muscle_HPA"), drop=FALSE]
colnames(DIA_Heart) <- c("GeneID","GeneName","Level")
colnames(HPA_Heart) <- c("GeneID","GeneName","Level")
Total_Heart <- unique(rbind(DIA_Heart[,c("GeneID"), drop=FALSE], HPA_Heart[,c("GeneID"), drop=FALSE]))
Common_Heart        <- subset(DIA_Heart, (GeneID %in% HPA_Heart$GeneID))
Only_in_DIA_Heart   <- subset(DIA_Heart, !(GeneID %in% HPA_Heart$GeneID))
Only_in_HPA_Heart   <- subset(HPA_Heart, !(GeneID %in% DIA_Heart$GeneID))

DIA_Liver <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Liver_DIA),c("GeneID_DIA","GeneName_DIA","Liver_DIA"), drop=FALSE]
HPA_Liver <- HPA_median_bins[!is.na(HPA_median_bins$liver_HPA),c("Gene_HPA","Gene.name_HPA","liver_HPA"), drop=FALSE]
colnames(DIA_Liver) <- c("GeneID","GeneName","Level")
colnames(HPA_Liver) <- c("GeneID","GeneName","Level")
Total_Liver <- unique(rbind(DIA_Liver[,c("GeneID"), drop=FALSE], HPA_Liver[,c("GeneID"), drop=FALSE]))
Common_Liver        <- subset(DIA_Liver, (GeneID %in% HPA_Liver$GeneID))
Only_in_DIA_Liver   <- subset(DIA_Liver, !(GeneID %in% HPA_Liver$GeneID))
Only_in_HPA_Liver   <- subset(HPA_Liver, !(GeneID %in% DIA_Liver$GeneID))

DIA_Lung <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Lung_DIA),c("GeneID_DIA","GeneName_DIA","Lung_DIA"), drop=FALSE]
HPA_Lung <- HPA_median_bins[!is.na(HPA_median_bins$lung_HPA),c("Gene_HPA","Gene.name_HPA","lung_HPA"), drop=FALSE]
colnames(DIA_Lung) <- c("GeneID","GeneName","Level")
colnames(HPA_Lung) <- c("GeneID","GeneName","Level")
Total_Lung <- unique(rbind(DIA_Lung[,c("GeneID"), drop=FALSE], HPA_Lung[,c("GeneID"), drop=FALSE]))
Common_Lung        <- subset(DIA_Lung, (GeneID %in% HPA_Lung$GeneID))
Only_in_DIA_Lung   <- subset(DIA_Lung, !(GeneID %in% HPA_Lung$GeneID))
Only_in_HPA_Lung   <- subset(HPA_Lung, !(GeneID %in% DIA_Lung$GeneID))

DIA_Pancreas <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Pancreas_DIA),c("GeneID_DIA","GeneName_DIA","Pancreas_DIA"), drop=FALSE]
HPA_Pancreas <- HPA_median_bins[!is.na(HPA_median_bins$pancreas_HPA),c("Gene_HPA","Gene.name_HPA","pancreas_HPA"), drop=FALSE]
colnames(DIA_Pancreas) <- c("GeneID","GeneName","Level")
colnames(HPA_Pancreas) <- c("GeneID","GeneName","Level")
Total_Pancreas <- unique(rbind(DIA_Pancreas[,c("GeneID"), drop=FALSE], HPA_Pancreas[,c("GeneID"), drop=FALSE]))
Common_Pancreas        <- subset(DIA_Pancreas, (GeneID %in% HPA_Pancreas$GeneID))
Only_in_DIA_Pancreas   <- subset(DIA_Pancreas, !(GeneID %in% HPA_Pancreas$GeneID))
Only_in_HPA_Pancreas   <- subset(HPA_Pancreas, !(GeneID %in% DIA_Pancreas$GeneID))

DIA_Thyroid <- DIA_3bins_gene[!is.na(DIA_3bins_gene$Thyroid_DIA),c("GeneID_DIA","GeneName_DIA","Thyroid_DIA"), drop=FALSE]
HPA_Thyroid <- HPA_median_bins[!is.na(HPA_median_bins$thyroid_gland_HPA),c("Gene_HPA","Gene.name_HPA","thyroid_gland_HPA"), drop=FALSE]
colnames(DIA_Thyroid) <- c("GeneID","GeneName","Level")
colnames(HPA_Thyroid) <- c("GeneID","GeneName","Level")
Total_Thyroid <- unique(rbind(DIA_Thyroid[,c("GeneID"), drop=FALSE], HPA_Thyroid[,c("GeneID"), drop=FALSE]))
Common_Thyroid        <- subset(DIA_Thyroid, (GeneID %in% HPA_Thyroid$GeneID))
Only_in_DIA_Thyroid   <- subset(DIA_Thyroid, !(GeneID %in% HPA_Thyroid$GeneID))
Only_in_HPA_Thyroid   <- subset(HPA_Thyroid, !(GeneID %in% DIA_Thyroid$GeneID))

# table with percentage genes common between MQ analysis and HPA
Percent_common <- data.frame("Percent" = (nrow(Common_Brain)/nrow(Total_Brain))*100, "Organ" = "Brain")
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Colon)/nrow(Total_Colon))*100, "Organ" = "Colon"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Duodenum)/nrow(Total_Duodenum))*100, "Organ" = "Duodenum"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Esophagus)/nrow(Total_Esophagus))*100, "Organ" = "Esophagus"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Heart)/nrow(Total_Heart))*100, "Organ" = "Heart"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Liver)/nrow(Total_Liver))*100, "Organ" = "Liver"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Lung)/nrow(Total_Lung))*100, "Organ" = "Lung"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Pancreas)/nrow(Total_Pancreas))*100, "Organ" = "Pancreas"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Thyroid)/nrow(Total_Thyroid))*100, "Organ" = "Thyroid"))
Percent_common$Present <- rep("Both in this study and HPA", nrow(Percent_common))

# table with percentage genes only identified in DIA analysis and NOT in HPA
Percent_only_in_DIA <- data.frame("Percent" = (nrow(Only_in_DIA_Brain)/nrow(Total_Brain))*100, "Organ" = "Brain")
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Colon)/nrow(Total_Colon))*100, "Organ" = "Colon"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Duodenum)/nrow(Total_Duodenum))*100, "Organ" = "Duodenum"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Esophagus)/nrow(Total_Esophagus))*100, "Organ" = "Esophagus"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Heart)/nrow(Total_Heart))*100, "Organ" = "Heart"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Liver)/nrow(Total_Liver))*100, "Organ" = "Liver"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Lung)/nrow(Total_Lung))*100, "Organ" = "Lung"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Pancreas)/nrow(Total_Pancreas))*100, "Organ" = "Pancreas"))
Percent_only_in_DIA <- rbind(Percent_only_in_DIA, data.frame("Percent" = (nrow(Only_in_DIA_Thyroid)/nrow(Total_Thyroid))*100, "Organ" = "Thyroid"))
Percent_only_in_DIA$Present <- rep("Present only in this study", nrow(Percent_only_in_DIA))

# table with percentage genes only identified in HPA and NOT in DIA analysis
Percent_only_in_HPA <- data.frame("Percent" = (nrow(Only_in_HPA_Brain)/nrow(Total_Brain))*100, "Organ" = "Brain")
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Colon)/nrow(Total_Colon))*100, "Organ" = "Colon"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Duodenum)/nrow(Total_Duodenum))*100, "Organ" = "Duodenum"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Esophagus)/nrow(Total_Esophagus))*100, "Organ" = "Esophagus"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Heart)/nrow(Total_Heart))*100, "Organ" = "Heart"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Liver)/nrow(Total_Liver))*100, "Organ" = "Liver"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Lung)/nrow(Total_Lung))*100, "Organ" = "Lung"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Pancreas)/nrow(Total_Pancreas))*100, "Organ" = "Pancreas"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Thyroid)/nrow(Total_Thyroid))*100, "Organ" = "Thyroid"))
Percent_only_in_HPA$Present <- rep("Present only in HPA", nrow(Percent_only_in_HPA))

plotdata <- rbind(Percent_common, Percent_only_in_DIA, Percent_only_in_HPA)

write.table(plotdata, "Comparison_of_genes_present_in_DIAanalysis_and_HPA.txt", sep = "\t", row.names = FALSE, quote = FALSE )
plotdata$Present <- gsub("Present only in", "Only in", plotdata$Present, perl=TRUE)


ggplot(plotdata, aes(x=Organ, y=Percent, fill=factor(Present,levels=c("Only in HPA","Only in this study","Both in this study and HPA")))) + 
  geom_bar(stat="identity") + 
  theme_bw()+
  labs(y="Percentage")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=12,),
        axis.title=element_text(size=12))+
  #scale_fill_manual(values = c("#440154FF","#2E6E8EFF","#FDE725FF"))  + 
  guides(fill=guide_legend(title="Present"))+
  ggtitle("Comparison of DIA data with genes identified from Human Protein Atlas (HPA)")

