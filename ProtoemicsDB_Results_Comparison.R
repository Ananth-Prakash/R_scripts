# To compare proteomics results with ProteomicsDB.

library(stringr)

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/ProteomicsDB")


ProtDB  <- read.table( "Proteomicsdb_API_response.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
ProtDB$normalised.intensity <- as.numeric(ProtDB$normalised.intensity)
ProtDB$normalised.intensity <- format(round(ProtDB$normalised.intensity, 3), nsmall = 2)
ProtDB$normalised.intensity <- as.numeric(ProtDB$normalised.intensity)
ProtDB <- ProtDB[,-c(1)]
ProtDB <- ProtDB[order(ProtDB$tissue),]

ProtDB <- unique(ProtDB)

ProtDB$Organ <- ProtDB$tissue
ProtDB$Organ <- gsub("cardia", "heart", ProtDB$Organ, perl=TRUE)
ProtDB$Organ <- gsub("cerebral cortex|prefrontal cortex", "brain", ProtDB$Organ, perl=TRUE)
ProtDB$Organ <- gsub("colon muscle", "colon", ProtDB$Organ, perl=TRUE)

ProtDB_filtered <- subset(ProtDB , Organ == "adrenal gland" | Organ == "adrenal gland"| Organ == "brain" | 
                             Organ == "colon" | Organ == "esophagus" | Organ == "gall bladder" | 
                             Organ == "heart" | Organ == "kidney" | Organ == "liver" | Organ == "lung"|
                             Organ == "lymph node" | Organ == "ovary" |Organ == "pancreas" | Organ == "placenta"|
                             Organ == "prostate gland" | Organ == "rectum" | Organ == "salivary gland" | 
                             Organ == "spleen" | Organ == "stomach" | Organ == "testis" | Organ == "thyroid gland"|
                             Organ == "tonsil" | Organ == "urinary bladder")

ProtDB_filtered <- ProtDB_filtered[,-c(2)]
ProtDB_filtered$Organ <- as.character(ProtDB_filtered$Organ)
ProtDB_filtered$UniProt <- as.character(ProtDB_filtered$UniProt)
#ProtDB_filtered$Organ <- factor(ProtDB_filtered$Organ)
ProtDB_filtered$Organ <- str_replace(ProtDB_filtered$Organ, "^\\w{1}", toupper)
ProtDB_filtered$Organ <- gsub("$","_ProteomicsDB", ProtDB_filtered$Organ, perl=TRUE)
ProtDB_filtered_aggregate <- aggregate(ProtDB_filtered[,2], by=list(ProtDB_filtered$UniProt, ProtDB_filtered$Organ), median, na.rm =TRUE)

colnames(ProtDB_filtered_aggregate) <- c("UniProt","Organ","normalised.intensity")
#
ProtDB_filtered_wide <- spread(ProtDB_filtered_aggregate, Organ, normalised.intensity)

#read median intensities from reanalysed datasets (ourstudy)
reanalysed_results <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-GeneNames-Median_intensities.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

reanalysed_results <- reanalysed_results[-c(34:37)]
colnames(reanalysed_results) <- gsub("_median", "_this study", colnames(reanalysed_results))

#read uniprot ENSG mappings
ID_maps <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/SupplementaryFile_3_UniProt_IDs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
ID_maps <- ID_maps[,c(1:3)]

reanalysed_results <- merge(x=ID_maps, y=reanalysed_results,
                            by.x=c("GeneID","GeneName"), by.y=c("GeneID","GeneName"), all.x=FALSE, all.y=FALSE)

Merged_data <- merge(x=reanalysed_results, y=ProtDB_filtered_wide,
                     by.x=c("UniProt_ID"), by.y=c("UniProt"), all.x=FALSE, all.y=FALSE)

AdrenalGland <- cbind(Merged_data[, c(1,2,3,grep("Adrenal", colnames(Merged_data)))], Organ="AdrenalGland")
colnames(AdrenalGland ) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Brain <- cbind(Merged_data[, c(1,2,3,grep("Brain", colnames(Merged_data)))], Organ="Brain")
colnames(Brain) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Colon <- cbind(Merged_data[, c(1,2,3,grep("Colon", colnames(Merged_data)))], Organ="Colon")
colnames(Colon) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Esophagus <- cbind(Merged_data[, c(1,2,3,grep("Esophagus", colnames(Merged_data)))], Organ="Esophagus")
colnames(Esophagus) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Gall_bladder <- cbind(Merged_data[, c(1,2,3,grep("Gall", colnames(Merged_data)))], Organ="Gallbladder")
colnames(Gall_bladder) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Heart <- cbind(Merged_data[, c(1,2,3,grep("Heart", colnames(Merged_data)))], Organ="Heart")
colnames(Heart) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Kidney <- cbind(Merged_data[, c(1,2,3,grep("Kidney", colnames(Merged_data)))], Organ="Kidney")
colnames(Kidney) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Liver <- cbind(Merged_data[, c(1,2,3,grep("Liver", colnames(Merged_data)))], Organ="Liver")
colnames(Liver) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Lung <- cbind(Merged_data[, c(1,2,3,grep("Lung", colnames(Merged_data)))], Organ="Lung")
colnames(Lung) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
LymphNode <- cbind(Merged_data[, c(1,2,3,grep("Lymph", colnames(Merged_data)))], Organ="LymphNode")
colnames(LymphNode) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Ovary <- cbind(Merged_data[, c(1,2,3,grep("Ovary", colnames(Merged_data)))], Organ="Ovary")
colnames(Ovary) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Pancreas <- cbind(Merged_data[, c(1,2,3,grep("Pancreas", colnames(Merged_data)))], Organ="Pancreas")
colnames(Pancreas) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Placenta <- cbind(Merged_data[, c(1,2,3,grep("Placenta", colnames(Merged_data)))], Organ="Placenta")
colnames(Placenta) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Prostate <- cbind(Merged_data[, c(1,2,3,grep("Prostate", colnames(Merged_data)))], Organ="Prostate")
colnames(Prostate) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Rectum <- cbind(Merged_data[, c(1,2,3,grep("Rectum", colnames(Merged_data)))], Organ="Rectum")
colnames(Rectum) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
SalivaryGland <- cbind(Merged_data[, c(1,2,3,grep("Salivary", colnames(Merged_data)))], Organ="SalivaryGland")
colnames(SalivaryGland) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Spleen <- cbind(Merged_data[, c(1,2,3,grep("Spleen", colnames(Merged_data)))], Organ="Spleen")
colnames(Spleen) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Stomach <- cbind(Merged_data[, c(1,2,3,grep("Stomach", colnames(Merged_data)))], Organ="Stomach")
colnames(Stomach) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Testis <- cbind(Merged_data[, c(1,2,3,grep("Testis", colnames(Merged_data)))], Organ="Testis")
colnames(Testis) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Thyroid <- cbind(Merged_data[, c(1,2,3,grep("Thyroid", colnames(Merged_data)))], Organ="Thyroid")
colnames(Thyroid) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
Tonsil <- cbind(Merged_data[, c(1,2,3,grep("Tonsil", colnames(Merged_data)))], Organ="Tonsil")
colnames(Tonsil) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")
UrinaryBladder <- cbind(Merged_data[, c(1,2,3,grep("Urinary", colnames(Merged_data)))], Organ="UrinaryBladder")
colnames(UrinaryBladder) <- c("UniProtID","GeneID","GeneName","Prakash_etal_2021", "ProteomicsDB", "Organ")



plotdata <- rbind(AdrenalGland, Brain, Colon, Esophagus, Gall_bladder, Heart, Liver, Lung, LymphNode, Ovary, Pancreas,
                  Placenta, Prostate, Rectum, SalivaryGland, Spleen, Stomach, Testis, Thyroid, Tonsil, UrinaryBladder)



ggplot(plotdata, aes(x=log2(Prakash_etal_2021), y=log2(ProteomicsDB))) + 
  geom_point(size=0.05, alpha=0.5) + 
  xlab("log2(FOT normalised iBAQ) this study")+
  ylab("log2(normalised intensity) ProteomicsDB")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  #scale_x_continuous(limits = c(0, 100000))+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Protein abundance comparison (This study vs ProteomicsDB)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  facet_wrap(~Organ)




