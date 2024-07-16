# To combine iBAQ Normal/Control samples on the Majority protein group ids
library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(mygene)

dataset1 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD012677/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset2 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD006692/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset3 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD016793/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset4 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD004364/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset5 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD001839/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset6 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD013543/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
dataset7 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD016958/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset8 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD003375/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)
dataset9 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD015928/proteinGroups-tissue_names-controlsonly.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)

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
  #tmp <- tmp[ ,-c(grep("iBAQ.peptides", colnames(tmp))) ]
  
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
dataset3 <- postprocess(dataset3)
dataset4 <- postprocess(dataset4)
dataset5 <- postprocess(dataset5)
dataset6 <- postprocess(dataset6)
dataset7 <- postprocess(dataset7)
dataset8 <- postprocess(dataset8)
dataset9 <- postprocess(dataset9)


#Remove trembl, swissprot (tr|sp) annotation before merging?
remove_ID_annotations <- function(dataframe){
  dataframe$Majority.protein.IDs <- gsub("sp\\||tr\\|","", dataframe$Majority.protein.IDs, perl=TRUE)
  dataframe$Majority.protein.IDs <- gsub("\\|\\w+_RAT","", dataframe$Majority.protein.IDs, perl=TRUE)
  return(dataframe)
}

dataset1 <- remove_ID_annotations(dataset1)
dataset2 <- remove_ID_annotations(dataset2)
dataset3 <- remove_ID_annotations(dataset3)
dataset4 <- remove_ID_annotations(dataset4)
dataset5 <- remove_ID_annotations(dataset5)
dataset6 <- remove_ID_annotations(dataset6)
dataset7 <- remove_ID_annotations(dataset7)
dataset8 <- remove_ID_annotations(dataset8)
dataset9 <- remove_ID_annotations(dataset9)

# Sort protein groups

dataset1$Majority.protein.IDs <- apply(dataset1, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset2$Majority.protein.IDs <- apply(dataset2, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset3$Majority.protein.IDs <- apply(dataset3, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset4$Majority.protein.IDs <- apply(dataset4, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset5$Majority.protein.IDs <- apply(dataset5, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset6$Majority.protein.IDs <- apply(dataset6, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset7$Majority.protein.IDs <- apply(dataset7, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset8$Majority.protein.IDs <- apply(dataset8, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))
dataset9$Majority.protein.IDs <- apply(dataset9, 1 , function(x) paste(sort(strsplit(x, ";")[[1]]), collapse = ";"))


Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Majority.protein.IDs"), by.y=c("Majority.protein.IDs"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(dataset1, dataset2)
merged_data <- Merge_data(merged_data, dataset3)
merged_data <- Merge_data(merged_data, dataset4)
merged_data <- Merge_data(merged_data, dataset5)
merged_data <- Merge_data(merged_data, dataset6)
merged_data <- Merge_data(merged_data, dataset7)
merged_data <- Merge_data(merged_data, dataset8)
merged_data <- Merge_data(merged_data, dataset9)




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
  res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("rat")),
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

merged_data <- merged_data[,c(1,99:102,2:98)]


merged_data_long <- gather(merged_data, Datasets, Intensities, colnames(merged_data[6]):colnames(merged_data[ncol(merged_data)]), factor_key=TRUE)

merged_data_long$Organs <- gsub("ppb.iBAQ.", "", merged_data_long$Datasets)

merged_data_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|Occipital|PinealGland|Pineal.Gland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex|MiddleBrain|PerfrontalCortex",
                                "Brain", merged_data_long$Organs, ignore.case=TRUE)
merged_data_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                "Heart", merged_data_long$Organs)
merged_data_long$Organs <- gsub("CaudalSegment|RostralSegment", "SpinalCord", merged_data_long$Organs)



merged_data_long$Organs <- gsub(".*Brain.*", "Brain", merged_data_long$Organs)

merged_data_long$Organs <- gsub("\\..*", "", merged_data_long$Organs)

merged_data_long_aggregate <- aggregate(merged_data_long[,-c(4,5,6)], list("MajorityProteinIDs" = merged_data_long$Majority.protein.IDs,
                                                                           "EnsemblID" = merged_data_long$ENSG,
                                                                           "GeneSymbol" = merged_data_long$Gene.Symbol,
                                                                           "Organs" = merged_data_long$Organs), median, na.rm =TRUE)

merged_data_long_aggregate <- merged_data_long_aggregate[,-c(5,6,7,9)]

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

Isoforms_organs_median_intensities <- Isoforms_organs_median_intensities[,c(1,2,4,3,13,5:12)]

Isoforms_organs_median_intensities <- Isoforms_organs_median_intensities[Isoforms_organs_median_intensities$present_in_samples > 0, ]

#To have the same number of genes in both files
median_bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Gene_distribution_in_organs-GeneNames-Median_bin_values-RAT.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)

Isoforms_organs_median_intensities <- merge(x=Isoforms_organs_median_intensities , y=median_bins[,c("GeneName"), drop=FALSE],
                                            by.x=c("Gene.Symbol"), by.y=c("GeneName"),
                                            all.x=FALSE, all.y=FALSE)

write.table(Isoforms_organs_median_intensities, file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Gene_distribution_in_organs-Isoforms-Median_intensities_new.txt", sep = "\t", row.names = FALSE, quote = FALSE )
