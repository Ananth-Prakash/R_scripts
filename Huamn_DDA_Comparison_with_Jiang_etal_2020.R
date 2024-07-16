#Comparison with Jiang 2020
# Supplementary table 2
library(tidyr)
library(ggpubr)
library(ggplot2)

Jiang_2020 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Normalised_protein_abundance_Jiang2020.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

colnames (Jiang_2020) <- gsub("\\.","", colnames(Jiang_2020), perl = TRUE)
#
Jiang_long <- gather(Jiang_2020, Organs, abundances, colnames(Jiang_2020)[3]:colnames(Jiang_2020)[ncol(Jiang_2020)], factor_key=TRUE)
Jiang_long$Organs <- gsub("\\d","", Jiang_long$Organs, perl = TRUE)
Jiang_long$Organs <- gsub("EsophagusGastroesophagealJunction|EsophagusMucosa|EsophagusMuscularis", "Esophagus", Jiang_long$Organs)
Jiang_long$Organs <- gsub("ArteryAorta|HeartAtrialAppendage|HeartLeftVentricle|ArteryCoronary", "Heart", Jiang_long$Organs)
Jiang_long$Organs <- gsub("ColonSigmoid|ColonTransverse", "Colon", Jiang_long$Organs)
Jiang_long$Organs <- gsub("BrainCerebellum|BrainCortex", "Brain", Jiang_long$Organs)
Jiang_long$Organs <- gsub("MinorSalivaryGland", "SalivaryGland", Jiang_long$Organs)
Jiang_long$Organs <- gsub("SmallIntestineTerminalIleum", "SmallIntestine", Jiang_long$Organs)

Jiang_long_aggregate <- aggregate(Jiang_long[,1:4], by=list(Jiang_long$geneidfull, Jiang_long$geneid, Jiang_long$Organs), median, na.rm =TRUE)

Jiang_long_aggregate <- Jiang_long_aggregate[,-c(4,5,6)]
colnames(Jiang_long_aggregate) <- c("GeneIDfull","GeneID","Organ","Abundance")
Jiang_long_aggregate <- Jiang_long_aggregate[!is.na(Jiang_long_aggregate$Abundance),]

Jiang_long_aggregate$Median_log2_Abundance <- log2(Jiang_long_aggregate$Abundance)

Jiang_long_aggregate <- Jiang_long_aggregate[,-c(4)]
Jiang_long_aggregate_wide <- spread(Jiang_long_aggregate, Organ, Median_log2_Abundance)
colnames(Jiang_long_aggregate_wide) <- gsub("$", "_Jiang_Etal_2020", colnames(Jiang_long_aggregate_wide), perl=TRUE)

#read median intensities from reanalysed datasets (ourstudy)
reanalysed_results <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-GeneNames-Median_intensities.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

reanalysed_results <- reanalysed_results[-c(34:37)]

reanalysed_results_long <- gather(reanalysed_results, Organs, Median_abundance, colnames(reanalysed_results)[3]:colnames(reanalysed_results)[ncol(reanalysed_results)], factor_key=TRUE)

reanalysed_results_long <- reanalysed_results_long[!is.na(reanalysed_results_long$Median_abundance),]
reanalysed_results_long$Median_log2_Abundance <- log2(reanalysed_results_long$Median_abundance)

reanalysed_results_long <- reanalysed_results_long[,-c(4)]

reanalysed_results_wide <- spread(reanalysed_results_long, Organs, Median_log2_Abundance)
colnames(reanalysed_results_wide) <- gsub("median", "Prakash_Etal_2021", colnames(reanalysed_results_wide), perl=TRUE)

Merged_data <- merge(x=reanalysed_results_wide, y=Jiang_long_aggregate_wide[-c(1)],
                     by.x=c("GeneID"), by.y=c("GeneID_Jiang_Etal_2020"), all.x=FALSE, all.y=FALSE)

AdrenalGland <- cbind(Merged_data[, c(1,2,grep("Adrenal", colnames(Merged_data)))], Organ="AdrenalGland")
colnames(AdrenalGland ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Brain <- cbind(Merged_data[, c(1,2,grep("Brain", colnames(Merged_data)))], Organ="Brain")
colnames(Brain ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Colon <- cbind(Merged_data[, c(1,2,grep("Colon", colnames(Merged_data)))], Organ="Colon")
colnames(Colon ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Esophagus <- cbind(Merged_data[, c(1,2,grep("Esophagus", colnames(Merged_data)))], Organ="Esophagus")
colnames(Esophagus ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Heart <- cbind(Merged_data[, c(1,2,grep("Heart", colnames(Merged_data)))], Organ="Heart")
colnames(Heart ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Liver <- cbind(Merged_data[, c(1,2,grep("Liver", colnames(Merged_data)))], Organ="Liver")
colnames(Liver ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Lung <- cbind(Merged_data[, c(1,2,grep("Lung", colnames(Merged_data)))], Organ="Lung")
colnames(Lung ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Ovary <- cbind(Merged_data[, c(1,2,grep("Ovary", colnames(Merged_data)))], Organ="Ovary")
colnames(Ovary ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Pancreas <- cbind(Merged_data[, c(1,2,grep("Pancreas", colnames(Merged_data)))], Organ="Pancreas")
colnames(Pancreas ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Prostate <- cbind(Merged_data[, c(1,2,grep("Prostate", colnames(Merged_data)))], Organ="Prostate")
colnames(Prostate ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
SalivaryGland <- cbind(Merged_data[, c(1,2,grep("SalivaryGland", colnames(Merged_data)))], Organ="SalivaryGland")
colnames(SalivaryGland ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
SmallIntestine <- cbind(Merged_data[, c(1,2,grep("SmallIntestine", colnames(Merged_data)))], Organ="SmallIntestine")
colnames(SmallIntestine ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Spleen <- cbind(Merged_data[, c(1,2,grep("Spleen", colnames(Merged_data)))], Organ="Spleen")
colnames(Spleen ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Stomach <- cbind(Merged_data[, c(1,2,grep("Stomach", colnames(Merged_data)))], Organ="Stomach")
colnames(Stomach ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Testis <- cbind(Merged_data[, c(1,2,grep("Testis", colnames(Merged_data)))], Organ="Testis")
colnames(Testis ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")
Thyroid <- cbind(Merged_data[, c(1,2,grep("Thyroid", colnames(Merged_data)))], Organ="Thyroid")
colnames(Thyroid ) <- c("GeneID","GeneName","Prakash_etal_2021", "Jiang_etal_2020", "Organ")

plotdata <- rbind(AdrenalGland, Brain, Colon, Esophagus, Heart, Liver, Lung, Ovary, Pancreas,
                  Prostate, SalivaryGland, SmallIntestine, Spleen, Stomach, Testis, Thyroid)

ggplot(plotdata, aes(x=Prakash_etal_2021, y=Jiang_etal_2020)) + 
  geom_point(size=0.05, alpha=0.5) + 
  xlab("log2(FOT normalised iBAQ) this study")+
  ylab("log2(normalised abundance) Jiang etal 2020")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Protein abundance comparison (This study vs Jiang et al 2020)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  facet_wrap(~Organ, scales="free_y")


