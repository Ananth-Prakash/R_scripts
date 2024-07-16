# To plot proteins identified vs/ MS spectra dectected (Quantify proteins vs. quantity of data) 
library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(ggpubr)

# Description of header
#http://www.coxdocs.org/doku.php?id=maxquant:table:summarytable

MS_spectra_data <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/All_DDA_datasets_Spectra_counts.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#MS_spectra_data <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/All_DDA_datasets_Spectra_counts_RAT.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)

MS_spectra_data$Dataset <- gsub(" ", "", MS_spectra_data$Dataset)

Canonical_protein_dist_datasets <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_datasets_plot.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#Canonical_protein_dist_datasets <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Gene_distribution_in_datasets_plot-RAT.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)


Canonical_protein_dist_organs <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs_plot.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#Canonical_protein_dist_organs <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Gene_distribution_in_organs_plot-RAT.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)


MS_spectra_data$Tissue <- gsub(".*_", "", MS_spectra_data$Experiment, perl=TRUE)

MS_spectra_data$Organs <- MS_spectra_data$Tissue
MS_spectra_data$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPituitaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|MiddleFrontalGyrus|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|Neocortex|OccipitalCortex|PinealGland|PituitaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                                "Brain", ignore.case = FALSE, MS_spectra_data$Organs)
MS_spectra_data$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                "Heart", ignore.case = FALSE, MS_spectra_data$Organs)

Organs_input <- subset(MS_spectra_data, select = c('MS', 'MS.MS', 'MS.MS.Identified', 'Organs'))
Spectra_sum_Organs <- aggregate(Organs_input[,c(1:3)], by=list(Organs=Organs_input$Organs), FUN=sum)
Spectra_sum_Organs <- Spectra_sum_Organs[Spectra_sum_Organs$Organs != "Breast",]
colnames(Spectra_sum_Organs) <- c("Organs", "number of MS spectra recorded", "number of MS/MS spectra recorded", "total number of identified\ntandem MS spectra")

Datasets_input <- subset(MS_spectra_data, select = c('MS', 'MS.MS', 'MS.MS.Identified', 'Dataset'))
Spectra_sum_Datasets <- aggregate(Datasets_input[,c(1:3)], by=list(Datasets=Datasets_input$Dataset), FUN=sum)
Spectra_sum_Datasets <- Spectra_sum_Datasets[Spectra_sum_Datasets$Datasets != "PXD001325",]
colnames(Spectra_sum_Datasets) <- c("Datasets", "number of MS spectra recorded", "number of MS/MS spectra recorded", "total number of identified\ntandem MS spectra")

Organs_MSspectra_CanonicalProtein_dist <- merge(x=Canonical_protein_dist_organs, y=Spectra_sum_Organs,
                                                by.x=c("Organs"), by.y=c("Organs"),all.x=TRUE, all.y=TRUE)

Datasets_MSspectra_CanonicalProtein_dist <- merge(x=Canonical_protein_dist_datasets, y=Spectra_sum_Datasets,
                                                by.x=c("Datasets"), by.y=c("Datasets"),all.x=TRUE, all.y=TRUE)

Organs_MSspectra_CanonicalProtein_dist_long <- gather(Organs_MSspectra_CanonicalProtein_dist, Index, Value, colnames(Organs_MSspectra_CanonicalProtein_dist)[6:8], factor_key=TRUE)
Organs_MSspectra_CanonicalProtein_dist_long$Type <- rep("Organs", nrow(Organs_MSspectra_CanonicalProtein_dist_long))

Datasets_MSspectra_CanonicalProtein_dist_long <- gather(Datasets_MSspectra_CanonicalProtein_dist, Index, Value, colnames(Datasets_MSspectra_CanonicalProtein_dist)[6:8], factor_key=TRUE)
Datasets_MSspectra_CanonicalProtein_dist_long$Type <- rep("Datasets", nrow(Datasets_MSspectra_CanonicalProtein_dist_long))


Organs_long_subset <- Organs_MSspectra_CanonicalProtein_dist_long[,c(2,5:8)]
colnames(Organs_long_subset)<- c("Number_of_identified_genes","Samples","Index","Value","Type")

Datasets_long_subset <- Datasets_MSspectra_CanonicalProtein_dist_long[,c(2,5:8)]
colnames(Datasets_long_subset)<- c("Number_of_identified_genes","Samples","Index","Value","Type")

plotdata <- rbind(Organs_long_subset, Datasets_long_subset)

#https://stackoverflow.com/questions/60143052/how-to-add-r2-for-each-facet-of-ggplot-in-r

#plotdata <- subset(plotdata, Index != "number of MS spectra recorded")

#plotdata <- plotdata[plotdata$Samples != "UmbilicalArtery (10)",]
#plotdata <- plotdata[plotdata$Samples != "PXD020187 (1)",]

ggplot(plotdata, aes(x=Value, y=Number_of_identified_genes)) + 
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
    facet_grid(Type~Index)
    #facet_grid(Type~Index, scales = "free_y")


