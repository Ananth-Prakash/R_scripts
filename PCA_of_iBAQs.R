#PCA of iBAQ values for response to reviewers
library(ggplot2)
library(dplyr)
library(stats)
library(tidyr)
library(viridis)
library(grid)
library(RColorBrewer)
library(heatmap3)
library(stringr)

#Brain
dataset1 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010154_Ananth/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD005819_33threads_yoda/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004143/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006233/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012755/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000547/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000548/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004332/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/ACT_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset10 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset11 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Banner_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset12 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_DorsoLateralPreFrontalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset13 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_Precuneus/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset14 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Mayo_TemporalCortex/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset15 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset16 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset17 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012131/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset18 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD015079/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset19 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010271/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset20 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006675/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset21 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008934/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")



Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene.ID", "Gene.Name"), by.y=c("Gene.ID", "Gene.Name"), all.x=TRUE, all.y=TRUE)
}

colnames(dataset1)[2] <- "Gene.Name"
colnames(dataset9)[2] <- "Gene.Name"
colnames(dataset10)[2] <- "Gene.Name"
colnames(dataset11)[2] <- "Gene.Name"
colnames(dataset12)[2] <- "Gene.Name"
colnames(dataset13)[2] <- "Gene.Name"
colnames(dataset14)[2] <- "Gene.Name"
colnames(dataset16)[2] <- "Gene.Name"
colnames(dataset17)[2] <- "Gene.Name"
colnames(dataset18)[2] <- "Gene.Name"


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

Merged_input_data[Merged_input_data == 0] <- NA


colnames(Merged_input_data) <- gsub("Gene.ID", "GeneID", colnames(Merged_input_data))
colnames(Merged_input_data) <- gsub("Gene.Name", "GeneName", colnames(Merged_input_data))
colnames(Merged_input_data) <- gsub(".*_", "", colnames(Merged_input_data), perl=TRUE)

Gene_info <- Merged_input_data[,c("GeneID", "GeneName")]
#Gene_info <- Gene_info[!duplicated(Gene_info$GeneID), ]

#Some of the gene entries (ex. IGHA2 has two Ensembl gene ids ENSG00000211890 & ENSG00000276173
#                          ex. IGHV2-70 has two Ensembl gene ids ENSG00000274576 & ENSG00000282453)
# because of this there are duplicate entries of such genes. These are aggregated by taking the median of them

Merged_input_data <- aggregate(Merged_input_data[,-c(1,2)], list("GeneName"= Merged_input_data$GeneName), median, na.rm =TRUE)

Brain <- Merged_input_data[grepl("GeneName|Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PituitaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex", colnames(Merged_input_data), ignore.case = TRUE)]
Heart <-  Merged_input_data[grepl("GeneName|Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum", colnames(Merged_input_data), ignore.case = TRUE)]


filtered_data <- Heart
#13168
# count number of samples for each gene that have non NA values (i.e., detected)
non_missing_sample_count_percentage <- data.frame(non_missing_sample_count_percentage= apply(filtered_data[2:ncol(filtered_data)], 1, function(x) (length(which(!is.na(x))))/(ncol(filtered_data)-1)*100 ))

filtered_data <- cbind(filtered_data, non_missing_sample_count_percentage)

# Filter genes that are present in at least 33%, 50% or 75% of samples
filtered_data <- filtered_data[filtered_data$non_missing_sample_count_percentage >= 50,]

### Since only common genes are considered, this condition is not required
# The undetected genes are also put into bin1 
# (After discussing with Andy Jones [20/01/2021]: since undetected genes could be below detection threshold due to their low expression or abundance)
filtered_data[is.na(filtered_data)] <- 1

#pca_input <- filtered_data[,-c(1, ncol(filtered_data), 30,31)]
pca_input <- filtered_data[,-c(1, ncol(filtered_data))]

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
                               var=unique(pca_plot_data$TissueID), 
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

ggplot(pca_plot_data, aes(x=PC1, y = PC2, colour = Tissues))+
  #geom_point(alpha=0.6)+
  geom_text(aes(label=TissueID))+
  labs(x="PC1", y="PC2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #labs(color="Samples") +
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.3,"line"))+
  guides(col = guide_legend(ncol = 3))+
  ggtitle("Brain-iBAQ\n[filter: genes detected in at least 50% of samples]\n[number of genes: 5357]")

# IMPORTANT reset key
GeomText$draw_key <- oldK
