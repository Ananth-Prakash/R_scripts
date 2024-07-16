#Compariaon of group and tissue enriched results with original published results PXD010154 (Supplementary table EV2)
# http://europepmc.org/article/MED/30777892#id246776
library(tidyverse)

PXD010154_data  <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/PXD010154_Supplementary_GroupTissue_Enriched.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

PXD010154_data$Tissue.enriched <- gsub("Fallopian tube", "FallopianTubeOviduct", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub(" gland", "Gland", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub(" intestine", "Intestine", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub("Urinary bladder", "UrinaryBladder", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub("Gallbladder", "GallBladder", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub("Appendix", "VermiformAppendix", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub("Fat", "AdiposeTissue", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub(" muscle", "Muscle", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub(" node", "Node", PXD010154_data$Tissue.enriched)
PXD010154_data$Tissue.enriched <- gsub("Endometrium", "UterineEndometrium", PXD010154_data$Tissue.enriched)


PXD010154_data$Group.enriched <- gsub("Fallopian tube", "FallopianTubeOviduct", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub(" gland", "Gland", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub(" intestine", "Intestine", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub("Urinary bladder", "UrinaryBladder", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub("Gallbladder", "GallBladder", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub("Appendix", "VermiformAppendix", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub("Fat", "AdiposeTissue", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub(" muscle", "Muscle", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub(" node", "Node", PXD010154_data$Group.enriched)
PXD010154_data$Group.enriched <- gsub("Endometrium", "UterineEndometrium", PXD010154_data$Group.enriched)



PXD010154_GroupEnriched <- PXD010154_data[PXD010154_data$Classification == "Group enriched", c("Gene.ID", "Gene.name", "Group.enriched", "Classification")]
colnames(PXD010154_GroupEnriched) <- c("Gene.ID", "Gene.name", "Organs", "Classification")
PXD010154_GroupEnriched$Organs <-  apply(PXD010154_GroupEnriched, 1 , function(x) paste(sort(strsplit(x, ";")[[3]]), collapse = ";"))


PXD010154_TissueEnriched <- PXD010154_data[PXD010154_data$Classification == "Tissue enriched", c("Gene.ID", "Gene.name", "Tissue.enriched", "Classification")]
colnames(PXD010154_TissueEnriched) <- c("Gene.ID", "Gene.name", "Organs", "Classification")


reanalysed_data <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/DDA_Manuscript_Supplementary_Files/oldfiles/SupplementaryFile_4.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
reanalysed_TissueEnriched <- reanalysed_data[reanalysed_data$Category == "Tissue enriched", ]


reanalysed_GroupEnriched <- reanalysed_data[reanalysed_data$Category == "Group enriched", ]


#reanalysed_data_grouped[reanalysed_data_grouped$Category=="Tissue enriched",]
#reanalysed_data_grouped <- reanalysed_data %>% group_by(Gene.ID, Gene.name, Category) %>%
#  summarise(Tissue = paste0(Tissue, collapse = ";"))

#reanalysed_data_grouped$Tissue <-  apply(reanalysed_data_grouped, 1 , function(x) paste(sort(strsplit(x, ";")[[4]]), collapse = ";"))
#colnames(reanalysed_data_grouped) <- c("Gene.ID","Gene.name","Classification","Organs")
#reanalysed_GroupEnriched <- reanalysed_data_grouped[reanalysed_data_grouped$Classification == "Group enriched", ]

# 1. Compare "Tissue enriched" set
merged_data_tissue_enriched <- merge(x=PXD010154_TissueEnriched, reanalysed_TissueEnriched,
                                     by.x=c("Gene.ID","Gene.name","Classification","Organs"),
                                     by.y=c("Gene.ID","Gene.name","Category","Tissue"),
                                     all.x=FALSE, all.y=FALSE)

tissue_enriched_count <- merged_data_tissue_enriched %>% count(Organs)


# 2. Compare "Group-enriched" set
PXD010154_Group_enriched_individuals <- PXD010154_GroupEnriched %>% 
  mutate(Organs=strsplit(Organs, ";")) %>% 
  unnest(Organs)

PXD010154_Group_enriched_individuals  <- unique(PXD010154_Group_enriched_individuals)

PXD010154_Group_enriched_wide <- spread(PXD010154_Group_enriched_individuals, Organs, Classification)
colnames(PXD010154_Group_enriched_wide) <- gsub("$","_PXD010154", colnames(PXD010154_Group_enriched_wide), perl=TRUE)

reanalysed_GroupEnriched_wide <- spread(reanalysed_GroupEnriched, Tissue, Category)
colnames(reanalysed_GroupEnriched_wide) <- gsub("$","_reanalysis", colnames(reanalysed_GroupEnriched_wide), perl=TRUE)

Merged_GroupEnriched <- merge(x=PXD010154_Group_enriched_wide, y=reanalysed_GroupEnriched_wide,
                              by.x=c("Gene.ID_PXD010154","Gene.name_PXD010154"),
                              by.y=c("Gene.ID_reanalysis","Gene.name_reanalysis"),
                              all.x=FALSE, all.y=FALSE)

group_enriched_common <- function(Organ){
    PXD_Tissue <- paste(Organ,"_PXD010154", sep="")
    reanalysis_Tissue <- paste(Organ,"_reanalysis", sep="")
    
    tmp <- Merged_GroupEnriched[,c("Gene.ID_PXD010154", "Gene.name_PXD010154", PXD_Tissue, reanalysis_Tissue)]
    
    tmp <- tmp[!is.na(tmp[[PXD_Tissue]]) | !is.na(tmp[[reanalysis_Tissue]]),]
    
    common <- tmp[!is.na(tmp[[PXD_Tissue]]) & !is.na(tmp[[reanalysis_Tissue]]),]
    percentage_common <- (nrow(common)/nrow(tmp))*100
    
    #print(head(common))
    print(percentage_common)
  return(percentage_common)
}


group_enrich_percent <- data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("AdiposeTissue"), Organ= "AdiposeTissue")
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("AdrenalGland"), Organ= "AdrenalGland"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Brain"), Organ= "Brain"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Colon"), Organ= "Colon"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Duodenum"), Organ= "Duodenum"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("UterineEndometrium"), Organ= "UterineEndometrium"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Esophagus"), Organ= "Esophagus"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("FallopianTubeOviduct"), Organ= "FallopianTubeOviduct"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("GallBladder"), Organ= "GallBladder"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Heart"), Organ= "Heart"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Kidney"), Organ= "Kidney"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Liver"), Organ= "Liver"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Lung"), Organ= "Lung"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("LymphNode"), Organ= "LymphNode"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Ovary"), Organ= "Ovary"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Pancreas"), Organ= "Pancreas"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Placenta"), Organ= "Placenta"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Prostate"), Organ= "Prostate"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Rectum"), Organ= "Rectum"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("SalivaryGland"), Organ= "SalivaryGland"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("SmallIntestine"), Organ= "SmallIntestine"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("SmoothMuscle"), Organ= "SmoothMuscle"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Spleen"), Organ= "Spleen"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Stomach"), Organ= "Stomach"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Testis"), Organ= "Testis"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Thyroid"), Organ= "Thyroid"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("Tonsil"), Organ= "Tonsil"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("UrinaryBladder"), Organ= "UrinaryBladder"))
group_enrich_percent <- rbind(group_enrich_percent, data.frame(Group_enriched_common_proteins_percentage=group_enriched_common("VermiformAppendix"), Organ= "VermiformAppendix"))

ggplot(group_enrich_percent, aes(x=Organ, y=Group_enriched_common_proteins_percentage)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("% of group enriched proteins common")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Group enriched proteins common between Prakash etal 2021 and PXD010154")


























