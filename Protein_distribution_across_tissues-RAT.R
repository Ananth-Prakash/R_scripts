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


#1. Separate samples by organs and extract genes
Brain <- Merged_input_data[grepl("GeneName|Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PituitaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex", colnames(Merged_input_data), ignore.case = TRUE)]
Heart <-  Merged_input_data[grepl("GeneName|Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum", colnames(Merged_input_data), ignore.case = TRUE)]
Kidney <- Merged_input_data[grepl("GeneName|Kidney|CorticalCollectingDuct|ConnectingTubule|CorticalThickAscendingLimb|DistalConvolutedTubule|InnerMedullaryCollectingDuct|MedullaryThickAscendingLimb|OuterMedullaryCollectingDuct|FirstSegmentOfProximalTubule|SecondSegmentOfProximalTubule|ThirdSegmentOfProximalTubule", colnames(Merged_input_data), ignore.case = TRUE)]
Liver <- Merged_input_data[grepl("GeneName|Liver", colnames(Merged_input_data), ignore.case = TRUE)]
Lung <- Merged_input_data[grepl("GeneName|Lung", colnames(Merged_input_data), ignore.case = TRUE)]
SpinalCord <- Merged_input_data[grepl("GeneName|SpinalCord|CaudalSegment|RostralSegment", colnames(Merged_input_data), ignore.case = TRUE)]
Tendon <- Merged_input_data[grepl("GeneName|Tendon", colnames(Merged_input_data), ignore.case = TRUE)]
Testis <- Merged_input_data[grepl("GeneName|Testis", colnames(Merged_input_data), ignore.case = TRUE)]

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



# Get sample sizes for each organ (i.e., number of MS runs representing each organ)
sample_sizes <- as.data.frame(ncol(Brain)-1)
sample_sizes <- cbind(sample_sizes, ncol(Heart)-1)
sample_sizes <- cbind(sample_sizes, ncol(Kidney)-1)
sample_sizes <- cbind(sample_sizes, ncol(Liver)-1)
sample_sizes <- cbind(sample_sizes, ncol(Lung)-1)
sample_sizes <- cbind(sample_sizes, ncol(SpinalCord)-1)
sample_sizes <- cbind(sample_sizes, ncol(Tendon)-1)
sample_sizes <- cbind(sample_sizes, ncol(Testis)-1)


colnames(sample_sizes) <- gsub(".*\\(", "", colnames(sample_sizes))
colnames(sample_sizes) <- gsub("\\).*", "", colnames(sample_sizes))
sample_sizes <- as.data.frame(t(sample_sizes))
sample_sizes <- tibble::rownames_to_column(sample_sizes, "Organs")
colnames(sample_sizes) <- c("Organs", "Number_of_samples")

# Data to plot distribution of iBAQ values across organs
Brain_iBAQ_long <- gather(Brain, Sample, iBAQ, colnames(Brain)[2]:colnames(Brain)[ncol(Brain)], factor_key=TRUE)
Heart_iBAQ_long <- gather(Heart, Sample, iBAQ, colnames(Heart)[2]:colnames(Heart)[ncol(Heart)], factor_key=TRUE)
Kidney_iBAQ_long <- gather(Kidney, Sample, iBAQ, colnames(Kidney)[2]:colnames(Kidney)[ncol(Kidney)], factor_key=TRUE)
Liver_iBAQ_long <- gather(Liver, Sample, iBAQ, colnames(Liver)[2]:colnames(Liver)[ncol(Liver)], factor_key=TRUE)
Lung_iBAQ_long <- gather(Lung, Sample, iBAQ, colnames(Lung)[2]:colnames(Lung)[ncol(Lung)], factor_key=TRUE)
Spinalcord_iBAQ_long <- gather(SpinalCord, Sample, iBAQ, colnames(SpinalCord)[2]:colnames(SpinalCord)[ncol(SpinalCord)], factor_key=TRUE)
Tendon_iBAQ_long <- gather(Tendon, Sample, iBAQ, colnames(Tendon)[2]:colnames(Tendon)[ncol(Tendon)], factor_key=TRUE)
Testis_iBAQ_long <- gather(Testis, Sample, iBAQ, colnames(Testis)[2]:colnames(Testis)[ncol(Testis)], factor_key=TRUE)


All_iBAQ_long <- rbind(Brain_iBAQ_long,Heart_iBAQ_long,Kidney_iBAQ_long,Liver_iBAQ_long,Lung_iBAQ_long,
                       Spinalcord_iBAQ_long,Tendon_iBAQ_long,Testis_iBAQ_long)

All_iBAQ_long <- All_iBAQ_long[complete.cases(All_iBAQ_long),]
All_iBAQ_long$Tissues <- gsub(".*\\.", "", All_iBAQ_long$Sample, perl=TRUE)

All_iBAQ_long$Organs <- All_iBAQ_long$Tissues
All_iBAQ_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPituitaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|MiddleFrontalGyrus|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|Neocortex|OccipitalCortex|PinealGland|PituitaryHypophysis|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                             "Brain", ignore.case = TRUE, All_iBAQ_long$Organs)
All_iBAQ_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                             "Heart", ignore.case = TRUE, All_iBAQ_long$Organs)
All_iBAQ_long$Organs  <- gsub("CaudalSegment1|CaudalSegment2|CaudalSegment3|RostralSegment1|RostralSegment2|RostralSegment3","SpinalCord", All_iBAQ_long$Organs , ignore.case=T, perl=TRUE)
All_iBAQ_long$Organs  <- gsub("Connectingtubule|Corticalcollectingduct|Corticalthickascendinglimb|Distalconvolutedtubule|Firstsegmentofproximaltubule|Innermedullarycollectingduct|Medullarythickascendinglimb|Outermedullarycollectingduct|Secondsegmentofproximaltubule|Thirdsegmentofproximaltubule","Kidney", All_iBAQ_long$Organs, ignore.case=T, perl=TRUE)


All_iBAQ_long <- merge(x=All_iBAQ_long, y=sample_sizes,
                       by.x=c("Organs"), by.y=c("Organs"))

All_iBAQ_long$Organs_samples <- paste(All_iBAQ_long$Organs, " (", All_iBAQ_long$Number_of_samples, ")", sep="")

All_iBAQ_long$Datasets <- gsub("\\..*", "", All_iBAQ_long$Sample, perl=TRUE)

ggplot(All_iBAQ_long, aes(x=Organs_samples, y=iBAQ)) + 
  geom_boxplot() + 
  xlab("Organs")+
  ylab("protein abundance (ppb)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("iBAQ (Rat)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))



# I. Compute the median from all samples of each Organ separately
Organs_median <- as.data.frame(cbind(GeneName = Brain$GeneName, Brain_median = apply(as.data.frame(Brain[,-c(1)]), 1, median, na.rm = T)))
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Heart$GeneName, Heart_median = apply(Heart[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Kidney$GeneName, Kidney_median = apply(as.data.frame(Kidney[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Liver$GeneName, Liver_median = apply(Liver[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Lung$GeneName, Lung_median = apply(as.data.frame(Lung[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = SpinalCord$GeneName, SpinalCord_median = apply(as.data.frame(SpinalCord[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Tendon$GeneName, Tendon_median = apply(as.data.frame(Tendon[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Organs_median <- merge(x=Organs_median, y=cbind(GeneName = Testis$GeneName, Testis_median = apply(as.data.frame(Testis[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)

write.table(Organs_median, file = paste("Gene_distribution_in_organs-GeneNames-Median_intensities-RAT.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

# II. Compute the median from all samples of each Tissue separately
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

Count <- data.frame(Present_in_number_of_samples = apply(Organs_median[2:ncol(Organs_median)], 1, function(x) length(which(x != "NaN")) ))
Count$Sample_percentage <- ((Count$Present_in_number_of_samples)/(ncol(Organs_median)-1)) *100
Count <- group_by(Count, Present_in_number_of_samples) %>% mutate(Total_number_of_sample_occurences = n()) %>% mutate(Gene_percent = (Total_number_of_sample_occurences / nrow(Count))*100)

Organs_median <- cbind(Organs_median, Count)

Organs_median <- merge(x=Gene_info, y=Organs_median,
                       by.x=c("GeneName"), by.y=c("GeneName"),
                       all.x=FALSE, all.y=FALSE)

plotdata <- unique(Count)


ggplot(plotdata[plotdata$Present_in_number_of_samples !=0,], aes(x=Present_in_number_of_samples, y=Gene_percent)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8"))+
  xlab("Present in number of organs")+
  ylab("% of identified 'canonical proteins'")+
  theme_bw()+
  ggtitle("Distribution of 'canonical proteins' across organs - (Rat)")

ggplot(plotdata[plotdata$Present_in_number_of_samples !=0,], aes(x=Present_in_number_of_samples, y=Total_number_of_sample_occurences)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8"))+
  xlab("Present in number of organs")+
  ylab("Number of canonical proteins")+
  theme_bw()+
  ggtitle("Distribution of 'canonical proteins' across organs (Rat)")


###### Gene distribution per organ
gene_counts_per_organ <- as.data.frame(sapply(Organs_median[,-c(1:2)], function(x) sum(!is.na(x))))
colnames(gene_counts_per_organ) <- "Number_of_identified_genes"
gene_counts_per_organ <- tibble::rownames_to_column(gene_counts_per_organ, "Organs")
gene_counts_per_organ$Organs <- gsub("_median", "", gene_counts_per_organ$Organs, perl=TRUE)
gene_counts_per_organ$Percentage <- (gene_counts_per_organ$Number_of_identified_genes/nrow(Organs_median))*100

gene_counts_per_organ <- merge(x=gene_counts_per_organ, y=sample_sizes,
                               by.x=c("Organs"), by.y="Organs")


gene_counts_per_organ$Organs_samples <- paste(gene_counts_per_organ$Organs, " (", gene_counts_per_organ$Number_of_samples, ")", sep="")

write.table(gene_counts_per_organ, file = paste("Gene_distribution_in_organs_plot-RAT.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

ggplot(gene_counts_per_organ, aes(x=Organs_samples, y=Percentage)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("% identified canonical proteins")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs (Rat)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggplot(gene_counts_per_organ, aes(x=Organs_samples, y=Number_of_identified_genes)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs (Rat)")+
  theme(axis.text=element_text(size=12,),
        axis.title=element_text(size=12))

##### Gene distribution per dataset
All_datasets <- Merged_input_data

dataset_sample_names <- colnames(All_datasets[-c(1)])
dataset_sample_names <- gsub("CaudalSegment1|CaudalSegment2|CaudalSegment3", "CaudalSegment", dataset_sample_names)
dataset_sample_names <- gsub("RostralSegment1|RostralSegment2|RostralSegment3", "RostralSegment", dataset_sample_names)
colnames(All_datasets) <- gsub("\\..*", "", colnames(All_datasets), perl=TRUE)


Datasets_median <- as.data.frame(cbind(GeneName = All_datasets$GeneName, PXD001839_median = apply(as.data.frame(All_datasets[,grep("PXD001839", colnames(All_datasets))]), 1, median, na.rm = T)))
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD003375_median = apply(as.data.frame(All_datasets[,grep("PXD003375", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD004364_median = apply(as.data.frame(All_datasets[,grep("PXD004364", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD006692_median = apply(as.data.frame(All_datasets[,grep("PXD006692", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD012677_median = apply(as.data.frame(All_datasets[,grep("PXD012677", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD013543_median = apply(as.data.frame(All_datasets[,grep("PXD013543", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD015928_median = apply(as.data.frame(All_datasets[,grep("PXD015928", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD016793_median = apply(as.data.frame(All_datasets[,grep("PXD016793", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD016958_median = apply(as.data.frame(All_datasets[,grep("PXD016958", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)


# Gene distribution per dataset
gene_counts_per_dataset <- as.data.frame(sapply(Datasets_median[,-c(1)], function(x) sum(!is.na(x))))
colnames(gene_counts_per_dataset) <- "Number_of_identified_genes"
gene_counts_per_dataset <- tibble::rownames_to_column(gene_counts_per_dataset, "Datasets")
gene_counts_per_dataset$Datasets <- gsub("_median", "", gene_counts_per_dataset$Datasets, perl=TRUE)
gene_counts_per_dataset$Percentage <- (gene_counts_per_dataset$Number_of_identified_genes/nrow(Datasets_median))*100

# Number of tissues per dataset
dataset_tissues <- data.frame(Datasets= "PXD001839", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD001839", dataset_sample_names)], perl=TRUE))) )
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD003375", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD003375", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD004364", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD004364", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD006692", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD006692", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD012677", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD012677", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD013543", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD013543", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD015928", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD015928", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD016793", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD016793", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD016958", Tissue_count = length(unique(gsub(".*\\.", "", dataset_sample_names[grep("PXD016958", dataset_sample_names)], perl=TRUE))) ))



gene_counts_per_dataset <- merge(x=gene_counts_per_dataset, y=dataset_tissues,
                                 by.x=c("Datasets"), by.y=c("Datasets"), all.x=FALSE, all.y=FALSE)

gene_counts_per_dataset$Datasets_tissues <- paste(gene_counts_per_dataset$Datasets, " (", gene_counts_per_dataset$Tissue_count, ")", sep="")

write.table(gene_counts_per_dataset, file = paste("Gene_distribution_in_datasets_plot-RAT.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )


ggplot(gene_counts_per_dataset, aes(x=Datasets_tissues, y=Number_of_identified_genes)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across datasets (Rat)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

# Plot distribution of iBAQ values across datasets
All_iBAQ_long <- merge(x=All_iBAQ_long, y=dataset_tissues,
                       by.x=c("Datasets"), by.y=c("Datasets"))

All_iBAQ_long$Dataset_tissues <- paste(All_iBAQ_long$Datasets, " (", All_iBAQ_long$Tissue_count, ")", sep="")


ggplot(All_iBAQ_long, aes(x=Dataset_tissues, y=iBAQ)) + 
  geom_boxplot() + 
  xlab("Datasets")+
  ylab("protein abundance (ppb)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("iBAQ (Rat)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))


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


###################################################################
############################## E N D ##############################
###################################################################
##### don't use the PCA here, use the method in rank-binning script  



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
