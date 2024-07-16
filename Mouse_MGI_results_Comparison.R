#Comparison with Mouse MGI
# http://www.informatics.jax.org/gxd

MGI <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/MGI/MGIgeneExpressionQuery_filtered.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)
MGI <- MGI[,c("Gene.Symbol","Structure","TPM.Level..RNA.Seq.")]

# To match bins
MGI$TPM.Level..RNA.Seq.[MGI$TPM.Level..RNA.Seq. == "Low"] <- 1
MGI$TPM.Level..RNA.Seq.[MGI$TPM.Level..RNA.Seq. == "Medium"] <- 2
MGI$TPM.Level..RNA.Seq.[MGI$TPM.Level..RNA.Seq. == "High"] <- 3

MGI$Gene.Symbol <- as.character(MGI$Gene.Symbol)
MGI$TPM.Level..RNA.Seq. <- as.numeric(MGI$TPM.Level..RNA.Seq.)

MGI_aggregate <- aggregate(MGI$TPM.Level..RNA.Seq., by = list(MGI$Gene.Symbol, MGI$Structure), FUN = median, na.rm =TRUE)
colnames(MGI_aggregate) <- c("GeneSymbol","Organ","Median_bins")

MGI_aggregate <- MGI_aggregate[MGI_aggregate$Organ != "brainstem",]
MGI_aggregate <- MGI_aggregate[MGI_aggregate$Organ != "heart left ventricle",]
MGI_aggregate <- MGI_aggregate[MGI_aggregate$Organ != "heart ventricle",]
MGI_aggregate <- MGI_aggregate[MGI_aggregate$Organ != "midbrain",]
MGI_aggregate <- MGI_aggregate[MGI_aggregate$Organ != "right lung middle lobe",]

MGI_aggregate$Organ <- gsub("brain","Brain_MGI", MGI_aggregate$Organ)
MGI_aggregate$Organ <- gsub("eye","Eye_MGI", MGI_aggregate$Organ)
MGI_aggregate$Organ <- gsub("heart","Heart_MGI", MGI_aggregate$Organ)
MGI_aggregate$Organ <- gsub("liver","Liver_MGI", MGI_aggregate$Organ)
MGI_aggregate$Organ <- gsub("lung","Lung_MGI", MGI_aggregate$Organ)
MGI_aggregate$Organ <- gsub("spleen","Spleen_MGI", MGI_aggregate$Organ)
MGI_aggregate$Organ <- gsub("testis","Testis_MGI", MGI_aggregate$Organ)

MGI_aggregate$Organ <- as.character(MGI_aggregate$Organ)
MGI_aggregate_wide <- spread(MGI_aggregate, Organ, Median_bins)

MQ_Mouse_3bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/Gene_distribution_in_organs-GeneNames-Median_3bins-Mouse.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(MQ_Mouse_3bins) <- gsub("_median", "_Wang.etal_2021", colnames(MQ_Mouse_3bins))

#MQ_Mouse_3bins_Genes <- MQ_Mouse_3bins[,c("GeneName","GeneID")]
#MQ_Mouse_3bins <- MQ_Mouse_3bins[,-c(1,2)]
#MQ_Mouse_3bins <- as.numeric(MQ_Mouse_3bins)

##################
#Spearman_correlation
##################
Spearman_correlation <- function(Tissue_MQ, Tissue_MGI){
  Merged_data <- merge(x = MQ_Mouse_3bins[,c("GeneName", Tissue_MQ)], 
                       y = MGI_aggregate_wide[,c("GeneSymbol", Tissue_MGI)],
                       by.x=c("GeneName"), by.y=c("GeneSymbol"), all.x=FALSE, all.y=FALSE)
  
  Merged_data[[Tissue_MQ]] <- as.numeric(Merged_data[[Tissue_MQ]])
  Merged_data[[Tissue_MGI]] <- as.numeric(Merged_data[[Tissue_MGI]])
  
  Merged_data <- Merged_data[!is.na(Merged_data[[Tissue_MQ]]) & !is.na(Merged_data[[Tissue_MGI]]),]
  
  #print(head(Merged_data))
  
  res <- cor.test(x=Merged_data[[Tissue_MQ]], y=Merged_data[[Tissue_MGI]], method = 'spearman')
  return(res$estimate[[1]])
} 

Spearman_correl <- data.frame(Spearmans_rho=Spearman_correlation("Brain_Wang.etal_2021", "Brain_MGI"), Organ="Brain")
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Eye_Wang.etal_2021", "Eye_MGI"), Organ="Eye"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Heart_Wang.etal_2021", "Heart_MGI"), Organ="Heart"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Liver_Wang.etal_2021", "Liver_MGI"), Organ="Liver"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Lung_Wang.etal_2021", "Lung_MGI"), Organ="Lung"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Spleen_Wang.etal_2021", "Spleen_MGI"), Organ="Spleen"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Testis_Wang.etal_2021", "Testis_MGI"), Organ="Testis"))

# MaxQuant+MGI vs MaxQuant+MGI
##############################
MQ_Mouse_3bins <- MQ_Mouse_3bins[-c(15:18)]
MQ_Mouse_3bins <- MQ_Mouse_3bins[MQ_Mouse_3bins$GeneName != "",]

combined_MaxQuant_MGI_bins <-merge(x = MQ_Mouse_3bins,y = MGI_aggregate_wide,
                                   by.x=c("GeneName"), by.y=c("GeneSymbol"), all.x=FALSE, all.y=FALSE)

edit_distances_combined_MaxQuant_MGI <- data.frame()

for(i in 3:ncol(combined_MaxQuant_MGI_bins)){
  for(j in 3:ncol(combined_MaxQuant_MGI_bins)){
    
    tmp <- data.frame()
    
    tmp <- merge(x = combined_MaxQuant_MGI_bins[,c("GeneName", colnames(combined_MaxQuant_MGI_bins)[i])] , 
                 y =  combined_MaxQuant_MGI_bins[,c("GeneName", colnames(combined_MaxQuant_MGI_bins)[j])] ,
                 by.x=c("GeneName"), by.y=c("GeneName"), all.x=FALSE, all.y=FALSE)
    
    tmp[[colnames(tmp)[2]]] <- as.numeric(tmp[[colnames(tmp)[2]]])
    tmp[[colnames(tmp)[3]]] <- as.numeric(tmp[[colnames(tmp)[3]]])
    
    tmp <- tmp[!is.na(tmp[[colnames(tmp)[2]]]) & !is.na(tmp[[colnames(tmp)[3]]]),]
    
    print(head(tmp))
    
    tmp$absdiff_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[3]]])
    average_editdist <- mean(tmp$absdiff_editdistance, na.rm=TRUE)
    
    result <- data.frame(Tissue1=colnames(tmp)[2], Tissue2=colnames(tmp)[3], average_absolute_editdistance=average_editdist)
    
    edit_distances_combined_MaxQuant_MGI <- rbind(edit_distances_combined_MaxQuant_MGI, result)
  }
}

edit_distances_combined_MaxQuant_MGI <- edit_distances_combined_MaxQuant_MGI[grep("Brain|Eye|Heart|Liver|Lung|Spleen|Testis", edit_distances_combined_MaxQuant_MGI$Tissue1, ignore.case = TRUE), ]
edit_distances_combined_MaxQuant_MGI <- edit_distances_combined_MaxQuant_MGI[grep("Brain|Eye|Heart|Liver|Lung|Spleen|Testis", edit_distances_combined_MaxQuant_MGI$Tissue2, ignore.case = TRUE), ]

edit_distances_combined_MaxQuant_MGI$Tissue1 <- gsub("\\.x","", edit_distances_combined_MaxQuant_MGI$Tissue1)
edit_distances_combined_MaxQuant_MGI$Tissue2 <- gsub("\\.y","", edit_distances_combined_MaxQuant_MGI$Tissue2)

edit_distances_combined_MaxQuant_MGI <- spread(edit_distances_combined_MaxQuant_MGI, Tissue2, average_absolute_editdistance)
rownames(edit_distances_combined_MaxQuant_MGI) <- edit_distances_combined_MaxQuant_MGI[,1]

edit_distances_combined_MaxQuant_MGI <- data.matrix(edit_distances_combined_MaxQuant_MGI[,-c(1)])

write.table(edit_distances_combined_MaxQuant_MGI, "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/MGI/Edit_distance_matrix_Mouse_Alldatasets_vs_Alldatasets.txt", sep = "\t", row.names = TRUE, quote = FALSE )

corrplot(edit_distances_combined_MaxQuant_MGI, is.corr = FALSE, type="full", tl.col="Black", tl.cex = 0.7, col = COL2('RdBu', 10))

par(oma=c(5,0,5,5)); 
heatmap.2(edit_distances_combined_MaxQuant_MGI, scale = "none", col = bluered(100), key.xlab="Edit distance",
          trace = "none", density.info = "none", cexRow=0.5, cexCol=0.75, margins = c(3,3))
