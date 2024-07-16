# To combine iBAQ Normal/Control samples on the Majority protein group ids
library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(mygene)

dataset1 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000547/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset2 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD000548/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#dataset3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001325/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset4 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001608_30threads_yoda/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD002029/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset6 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004143/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset7 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004332/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset8 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD005819_33threads_yoda/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset9 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006233/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset10 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD006675/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset11 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008934/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset12 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010154_Ananth/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset13 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010271/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset14 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012131/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset15 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012755/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset16 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD015079/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset17 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD020187/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset18 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/ACT_DorsoLateralPreFrontalCortex/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset19 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset20 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Banner_DorsoLateralPreFrontalCortex/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset21 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_DorsoLateralPreFrontalCortex/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset22 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/BLSA_Precuneus/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset23 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/Mayo_TemporalCortex/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset24 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/MountSinai_FrontalPole/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset25 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#dataset26 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012431/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#dataset27 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008722/proteinGroups-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)

##############
# define FOT normalisation function
# FOT stands for Fraction Of Total. In this normalisation method each protein iBAQ intensity value is scaled to the total amount of signal in a given MS run (column) and transformed to parts per billion (ppb)
fot.normalise <- function(x){
  data.sum <-   apply(x, 2, function(y){sum(y, na.rm=TRUE)})
  # barplot((data.sum), log = "y")
  #
  ##### do ppm normalisation
  x.mat <- as.matrix(x)
  x.mat.ppb <- apply(x.mat, 2, function(i) i/sum(i, na.rm = T) * 1000000000 )
  x.mat.ppb <- as.data.frame(x.mat.ppb)
  colnames(x.mat.ppb) <- paste("ppb.", colnames(x.mat.ppb), sep = "")
  return(x.mat.ppb)
}
##############


postprocess <- function(tmp){
  
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]

## consider those protein groups with more than 1 peptides matched
tmp <- tmp[ tmp[, "Peptides"] > 1, ]

tmp <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]


Majority.protein.IDs <- tmp$Majority.protein.IDs
tmp <- tmp[ , -1]
tmp <- fot.normalise(tmp)  

tmp <- data.frame( cbind(Majority.protein.IDs, tmp, stringsAsFactors = FALSE) )
##############
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA

return(tmp)
}

dataset1 <- postprocess(dataset1)
dataset2 <- postprocess(dataset2)
#dataset3 <- postprocess(dataset3)
dataset4 <- postprocess(dataset4)
dataset5 <- postprocess(dataset5)
dataset6 <- postprocess(dataset6)
dataset7 <- postprocess(dataset7)
dataset8 <- postprocess(dataset8)
dataset9 <- postprocess(dataset9)
dataset10 <- postprocess(dataset10)
dataset11 <- postprocess(dataset11)
dataset12 <- postprocess(dataset12)
dataset13 <- postprocess(dataset13)
dataset14 <- postprocess(dataset14)
dataset15 <- postprocess(dataset15)
dataset16 <- postprocess(dataset16)
dataset17 <- postprocess(dataset17)
dataset18 <- postprocess(dataset18)
dataset19 <- postprocess(dataset19)
dataset20 <- postprocess(dataset20)
dataset21 <- postprocess(dataset21)
dataset22 <- postprocess(dataset22)
dataset23 <- postprocess(dataset23)
dataset24 <- postprocess(dataset24)
dataset25 <- postprocess(dataset25)
#dataset26 <- postprocess(dataset26)
#dataset27 <- postprocess(dataset27)

#Remove trembl, swissprot (tr|sp) annotation before merging?
remove_ID_annotations <- function(dataframe){
  dataframe$Majority.protein.IDs <- gsub("sp\\||tr\\|","", dataframe$Majority.protein.IDs, perl=TRUE)
  dataframe$Majority.protein.IDs <- gsub("\\|\\w+_HUMAN","", dataframe$Majority.protein.IDs, perl=TRUE)
  return(dataframe)
}

dataset1 <- remove_ID_annotations(dataset1)
dataset2 <- remove_ID_annotations(dataset2)
#dataset3 <- remove_ID_annotations(dataset3)
dataset4 <- remove_ID_annotations(dataset4)
dataset5 <- remove_ID_annotations(dataset5)
dataset6 <- remove_ID_annotations(dataset6)
dataset7 <- remove_ID_annotations(dataset7)
dataset8 <- remove_ID_annotations(dataset8)
dataset9 <- remove_ID_annotations(dataset9)
dataset10 <- remove_ID_annotations(dataset10)
dataset11 <- remove_ID_annotations(dataset11)
dataset12 <- remove_ID_annotations(dataset12)
dataset13 <- remove_ID_annotations(dataset13)
dataset14 <- remove_ID_annotations(dataset14)
dataset15 <- remove_ID_annotations(dataset15)
dataset16 <- remove_ID_annotations(dataset16)
dataset17 <- remove_ID_annotations(dataset17)
dataset18 <- remove_ID_annotations(dataset18)
dataset19 <- remove_ID_annotations(dataset19)
dataset20 <- remove_ID_annotations(dataset20)
dataset21 <- remove_ID_annotations(dataset21)
dataset22 <- remove_ID_annotations(dataset22)
dataset23 <- remove_ID_annotations(dataset23)
dataset24 <- remove_ID_annotations(dataset24)
dataset25 <- remove_ID_annotations(dataset25)
#dataset26 <- remove_ID_annotations(dataset26)
#dataset27 <- remove_ID_annotations(dataset27)

# Sort protein groups

dataset1$Majority.protein.IDs <- apply(dataset1, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset2$Majority.protein.IDs <- apply(dataset2, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
#dataset3$Majority.protein.IDs <- apply(dataset3, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset4$Majority.protein.IDs <- apply(dataset4, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset5$Majority.protein.IDs <- apply(dataset5, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset6$Majority.protein.IDs <- apply(dataset6, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset7$Majority.protein.IDs <- apply(dataset7, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset8$Majority.protein.IDs <- apply(dataset8, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset9$Majority.protein.IDs <- apply(dataset9, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset10$Majority.protein.IDs <- apply(dataset10, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset11$Majority.protein.IDs <- apply(dataset11, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset12$Majority.protein.IDs <- apply(dataset12, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset13$Majority.protein.IDs <- apply(dataset13, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset14$Majority.protein.IDs <- apply(dataset14, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset15$Majority.protein.IDs <- apply(dataset15, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset16$Majority.protein.IDs <- apply(dataset16, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset17$Majority.protein.IDs <- apply(dataset17, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset18$Majority.protein.IDs <- apply(dataset18, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset19$Majority.protein.IDs <- apply(dataset19, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset20$Majority.protein.IDs <- apply(dataset20, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset21$Majority.protein.IDs <- apply(dataset21, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset22$Majority.protein.IDs <- apply(dataset22, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset23$Majority.protein.IDs <- apply(dataset23, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset24$Majority.protein.IDs <- apply(dataset24, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset25$Majority.protein.IDs <- apply(dataset25, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
#dataset26$Majority.protein.IDs <- apply(dataset26, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
#dataset27$Majority.protein.IDs <- apply(dataset27, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))



Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Majority.protein.IDs"), by.y=c("Majority.protein.IDs"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(dataset1, dataset2)
#merged_data <- Merge_data(merged_data, dataset3)
merged_data <- Merge_data(merged_data, dataset4)
merged_data <- Merge_data(merged_data, dataset5)
merged_data <- Merge_data(merged_data, dataset6)
merged_data <- Merge_data(merged_data, dataset7)
merged_data <- Merge_data(merged_data, dataset8)
merged_data <- Merge_data(merged_data, dataset9)
merged_data <- Merge_data(merged_data, dataset10)
merged_data <- Merge_data(merged_data, dataset11)
merged_data <- Merge_data(merged_data, dataset12)
merged_data <- Merge_data(merged_data, dataset13)
merged_data <- Merge_data(merged_data, dataset14)
merged_data <- Merge_data(merged_data, dataset15)
merged_data <- Merge_data(merged_data, dataset16)
merged_data <- Merge_data(merged_data, dataset17)
merged_data <- Merge_data(merged_data, dataset18)
merged_data <- Merge_data(merged_data, dataset19)
merged_data <- Merge_data(merged_data, dataset20)
merged_data <- Merge_data(merged_data, dataset21)
merged_data <- Merge_data(merged_data, dataset22)
merged_data <- Merge_data(merged_data, dataset23)
merged_data <- Merge_data(merged_data, dataset24)
merged_data <- Merge_data(merged_data, dataset25)
#merged_data <- Merge_data(merged_data, dataset26)
#merged_data <- Merge_data(merged_data, dataset27)



merged_data$"ENSG" <- "NA"
merged_data$"Gene.Symbol" <- "NA"
merged_data$"unique.gene.count" <- "NA"
merged_data$"protein.count" <- "NA"

for(i in 1:nrow(merged_data)){
  
  x <- data.frame(strsplit(merged_data[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
  x_temp <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  #x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
  #x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  #x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
  f = file()
  sink(file=f)
  res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("human")),
                  error = function(e) {print(0)})
  sink()
  close(f)
  
  
  if (class(res)=="DFrame" | class(res) == "DataFrame"){
    merged_data[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    merged_data[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")}
  temp_symb <- merged_data[i,"Gene.Symbol"]
  merged_data[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
  merged_data[ i , "protein.count"] <- str_count(merged_data[i, "Majority.protein.IDs"], ";")+1
  
  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(merged_data),1)),"%"))
}

backup <- merged_data

merged_data <- merged_data[ merged_data$Gene.Symbol != "NA" , ]
merged_data <- merged_data[ merged_data[, "unique.gene.count"] == 1, ]

merged_data <- merged_data[,-c(88)]
merged_data <- merged_data[,c(1,501:504,2:500)]

#merged_data <- merged_data[,c(1,505:508,2:504)]
#merged_data <- merged_data[,-c(95)]

write.table(merged_data, file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/protein-groups-isoforms_backup_1.txt", sep = "\t", row.names = FALSE, quote = FALSE )

#merged_data <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/protein-groups-isoforms_backup_1.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)


merged_data_long <- gather(merged_data, Datasets, Intensities, colnames(merged_data[6]):colnames(merged_data[ncol(merged_data)]), factor_key=TRUE)

merged_data_long$Organs <- gsub("ppb.iBAQ.", "", merged_data_long$Datasets)

merged_data_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|Occipital|PinealGland|Pineal.Gland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex|MiddleBrain|PerfrontalCortex",
                                "Brain", merged_data_long$Organs, ignore.case=TRUE)
merged_data_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                "Heart", merged_data_long$Organs)
merged_data_long$Organs <- gsub(".*Brain.*", "Brain", merged_data_long$Organs)
merged_data_long$Organs <- gsub(".*Breast.*", "Breast", merged_data_long$Organs)
merged_data_long$Organs <- gsub("\\..*", "", merged_data_long$Organs)


#merged_data_long_aggregate <- aggregate(merged_data_long[,-c(4,5,6)], by=list("MajorityProteinIDs" = merged_data_long$Majority.protein.IDs,
#                                                                         "EnsemblID" = merged_data_long$ENSG,
#                                                                         "GeneSymbol" = merged_data_long$Gene.Symbol,
#                                                                         "Organs" = merged_data_long$Organs), median, na.rm =TRUE)

merged_data_long <- merged_data_long[,-c(4,5,6)]
merged_data_long_aggregate <- aggregate(merged_data_long$Intensities, by=list(merged_data_long$Majority.protein.IDs,
                                                       merged_data_long$ENSG,
                                                       merged_data_long$Gene.Symbol,
                                                       merged_data_long$Organs), median, na.rm =TRUE)

colnames(merged_data_long_aggregate) <- c("MajorityProteinIDs","EnsemblID","GeneSymbol","Organs","Intensities")

Isoforms_organs_median_intensities <- spread(merged_data_long_aggregate, Organs, Intensities)


# Look for proteins that occurs uniquely in a protein group and not present in any other protein groups
single_protein_group <- merged_data[merged_data$protein.count == 1 ,c("Majority.protein.IDs", "Gene.Symbol")]
single_protein_group$unique_isoform_protein_group <- NA

multi_protein_group <- merged_data [ !merged_data$Majority.protein.IDs %in% single_protein_group$Majority.protein.IDs ,]

for(i in 1:nrow(single_protein_group)){
  x_protein.ID <- single_protein_group[ i, "Majority.protein.IDs"]
  x_gene.symbol <- single_protein_group[ i, "Gene.Symbol"]
  
  sub_multi_protein_group <- multi_protein_group[multi_protein_group$Gene.Symbol == x_gene.symbol, c("Majority.protein.IDs"), drop=FALSE]
   
  if (dim(sub_multi_protein_group)[1] != 0) {
    
    y_protein.IDs <- strsplit(paste((pull(sub_multi_protein_group, Majority.protein.IDs)),collapse=";"), ";")
    
    if( all(is.na(str_match(x_protein.ID, y_protein.IDs[[1]]))) == "TRUE"){
        single_protein_group[ i, "unique_isoform_protein_group"] <- "TRUE"
        }# end of if
    else{single_protein_group[ i, "unique_isoform_protein_group"] <- "FALSE"}
  }# end of if
  
  else{
    single_protein_group[ i, "unique_isoform_protein_group"] <- "TRUE"
  }# end of else
    
  } # end of for

Isoforms_organs_median_intensities <- merge(x=single_protein_group, y=Isoforms_organs_median_intensities,
                                            by.x=c("Majority.protein.IDs", "Gene.Symbol"), by.y=c("MajorityProteinIDs", "GeneSymbol"),
                                            all.x=TRUE, all.y=TRUE)

Isoforms_organs_median_intensities <- Isoforms_organs_median_intensities %>% replace_na(list(unique_isoform_protein_group = "FALSE"))

Isoforms_organs_median_intensities$present_in_samples <- rowSums(!is.na(Isoforms_organs_median_intensities[, 5:ncol(Isoforms_organs_median_intensities)]))

Isoforms_organs_median_intensities <- Isoforms_organs_median_intensities[,c(1,2,4,3,36,5:35)]

Isoforms_organs_median_intensities <- Isoforms_organs_median_intensities[Isoforms_organs_median_intensities$present_in_samples > 0, ]

#To have the same number of genes in both files
median_bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-GeneNames-Median_bin_values.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
#isoforms <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-Isoforms-Median_intensities_new.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)


Isoforms_organs_median_intensities <- merge(x=Isoforms_organs_median_intensities , y=median_bins[,c("GeneName"), drop=FALSE],
                  by.x=c("Gene.Symbol"), by.y=c("GeneName"),
                  all.x=FALSE, all.y=FALSE)
  
write.table(Isoforms_organs_median_intensities, file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-Isoforms-Median_intensities_new.txt", sep = "\t", row.names = FALSE, quote = FALSE )



###### Get the numbers for publication
### Use this function for post-process
# postprocess <- function(tmp){

#  tmp <- tmp[ tmp[, "Reverse"] != "+", ]
#  tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]
#  return(tmp)
# }

Total_protein_groups <- nrow(dataset1)+nrow(dataset2)+nrow(dataset3)+nrow(dataset4)+nrow(dataset5)+nrow(dataset6)+nrow(dataset7)+nrow(dataset8)+nrow(dataset9)+nrow(dataset10)+nrow(dataset11)+nrow(dataset12)+nrow(dataset13)+nrow(dataset14)+nrow(dataset15)+nrow(dataset16)+nrow(dataset17)+nrow(dataset18)+nrow(dataset19)+nrow(dataset20)+nrow(dataset21)+nrow(dataset22)+nrow(dataset23)+nrow(dataset24)+nrow(dataset25)
#134740

Total_peptides <- sum(dataset1$Peptides)+sum(dataset2$Peptides)+sum(dataset3$Peptides)+sum(dataset4$Peptides)+sum(dataset5$Peptides)+sum(dataset6$Peptides)+sum(dataset7$Peptides)+sum(dataset8$Peptides)+sum(dataset9$Peptides)+sum(dataset10$Peptides)+sum(dataset11$Peptides)+sum(dataset12$Peptides)+sum(dataset13$Peptides)+sum(dataset14$Peptides)+sum(dataset15$Peptides)+sum(dataset16$Peptides)+sum(dataset17$Peptides)+sum(dataset18$Peptides)+sum(dataset19$Peptides)+sum(dataset20$Peptides)+sum(dataset21$Peptides)+sum(dataset22$Peptides)+sum(dataset23$Peptides)+sum(dataset24$Peptides)+sum(dataset25$Peptides)
#2478836

Total_unique_peptides <-  sum(dataset1$Unique.peptides)+sum(dataset2$Unique.peptides)+sum(dataset3$Unique.peptides)+sum(dataset4$Unique.peptides)+sum(dataset5$Unique.peptides)+sum(dataset6$Unique.peptides)+sum(dataset7$Unique.peptides)+sum(dataset8$Unique.peptides)+sum(dataset9$Unique.peptides)+sum(dataset10$Unique.peptides)+sum(dataset11$Unique.peptides)+sum(dataset12$Unique.peptides)+sum(dataset13$Unique.peptides)+sum(dataset14$Unique.peptides)+sum(dataset15$Unique.peptides)+sum(dataset16$Unique.peptides)+sum(dataset17$Unique.peptides)+sum(dataset18$Unique.peptides)+sum(dataset19$Unique.peptides)+sum(dataset20$Unique.peptides)+sum(dataset21$Unique.peptides)+sum(dataset22$Unique.peptides)+sum(dataset23$Unique.peptides)+sum(dataset24$Unique.peptides)+sum(dataset25$Unique.peptides)
#1428859

library(plyr)
moo <- plyr::rbind.fill(dataset1[,c("Majority.protein.IDs"), drop=FALSE], 
                        dataset2[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset3[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset4[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset5[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset6[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset7[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset8[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset9[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset10[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset11[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset12[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset13[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset14[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset15[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset16[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset17[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset18[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset19[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset20[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset21[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset22[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset23[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset24[,c("Majority.protein.IDs"), drop=FALSE],
                        dataset25[,c("Majority.protein.IDs"), drop=FALSE])

Total_unique_genes <- length(unique(sort(unlist(apply(moo, 1 , function(x) strsplit(x, ";")[[1]]))))) 
#93449