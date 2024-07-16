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

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas")


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

write.table(HPA_median_bins, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/HPA_3bins.txt", sep = "\t", row.names = FALSE, quote = FALSE )


#Compare genes identified in MQ analysis with data from HPA

#MQ_analysis  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-GeneNames-Median_bin_values.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#colnames(MQ_analysis) <- gsub("$", "_Manuscript_analysis", colnames(MQ_analysis), perl=TRUE)

#Andy's suggestion: Rebin with only 3 bins
####### NOT E#########################################################################
# Don't use this function, gives wrong binning results for those with only 1 sample
# Use excel rebinned data insted, see random edit distance at the bottom (line 1116)
####### NOT E#########################################################################
Three_bins <- function(datasetID, tissue){

  
  tmp  <- read.table(file = paste("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/", datasetID, "/proteinGroups_ppb_final-tissue_names.txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
  
  colnames(tmp)[2] <- "Gene.Symbol"
  dataset <- aggregate(tmp[,-c(1,2)], list(tmp$Gene.Symbol), median, na.rm =TRUE)
  colnames(dataset)[1] <- "Gene" 
  dataset[dataset==0] <- NA
  
  colnames(dataset) <- gsub(".*_", "", colnames(dataset), perl=TRUE)
  gene_names <- dataset[,c(1)] 
  
  
  convert_to_scaled_data <- function(data, bins) {
    apply(data, 1, function(sample) {
      ceiling((sample / max(sample, na.rm = T)) * bins)
    })
  }
  
  ######Original
  convert_to_ranked_data <- function(data) {
    apply(data, 2, function(sample_data){
      unique_vals = unique(sample_data)
      unique_vals = unique_vals[order(unique_vals)]
      sapply(sample_data, function(s){which(s == unique_vals)}) - 1
    })
  }
  
  #convert_to_ranked_data <- function(data) {
  #  apply(data, 2, function(sample_data){
  #    unique_vals = unique(sample_data)
  #    #unique_vals = unique_vals[order(unique_vals)]
  #    unique_vals = unique_vals[order(unique_vals, na.last = NA)]
  #    #sapply(sample_data, function(s){which(s == unique_vals)}) - 1
  #    sapply(sample_data, function(s){ if(!is.na(s)){which(s == unique_vals)} else {s=NA} }) - 1
  #  })
  #}
  
  Binning <- function(data_input, gene_names_inp) {
    # Make a data frame of sample + protein/peptide abundances. 
    data = data.frame(data_input)
    
    if( ncol(data) == 1){
      data <- cbind(data, gene_names_inp) 
      data <- data.frame(data[complete.cases(data),])
      gene_names <- data[,c("gene_names_inp")]
      data <- subset(data, select=-c(gene_names_inp))
    }
    
    # Rank data.
    data_ranked = as.data.frame(t(convert_to_ranked_data(data)))
    colnames(data_ranked) <- gene_names
    
    # Different size of bins to try.
    binning_size_approaches = list(
      bins_3 = 2
      #bins_5 = 4,
      #bins_10 = 9,
      #bins_100 = 99,
      #bins_1000 = 999,
      #bins_10000 = 9999
    )
    
    # Transform data using different bin sizes.
    results = lapply(binning_size_approaches, function(binning_size_approach) {
      # 'Bin' data by scaling the ranked data.
      binned_data = as.data.frame(t(convert_to_scaled_data(data_ranked, binning_size_approach)))
    })
  }
  
  # If a dataset has different tissues or tissue regions, do not bin all together,
  # these tissues/regions must be separated and binned separately.
  #binning_input <- dataset[,-c(1)]
  binning_input <- dataset[,c(grep(tissue, colnames(dataset))), drop=FALSE]
  
  #keep a copy of the input before changing NAs to 0
  binning_input_with_NA <- binning_input
  #binning_input_with_NA <- dataset
  
  if(ncol(binning_input) > 1){
    # non-detected genes (NA or 0 values) are assigned 0 values here to help in binning.
    binning_input[is.na(binning_input)] <- 0
  }
  
  Dataset_all_bins <- Binning(binning_input, gene_names)
  
  # Consider only data grouped into 3 bins
  Get_Bin_3 <- function(binned_data_list){
    
    Bins_3 <- as.data.frame(t(binned_data_list$bins_3))
    # rename bins 0 to 2 --> 1 to 3
    # lowest intensities are in bin1 and highest in bin3
    Bins_3[Bins_3 == 2] <- 3
    Bins_3[Bins_3 == 1] <- 2
    Bins_3[Bins_3 == 0] <- 1
    
    Bins_3<- tibble::rownames_to_column(Bins_3, "Gene")
  }
  
  Dataset_3_bins <- Get_Bin_3(Dataset_all_bins)
  
  # Now reassign NA values back to those genes which were undetected or 0 in ALL samples 
  # NOTE: Only those genes are reassigned back to NA that have missing values in
  # all samples (ie., if at least one sample had a value then it is left as 0)

  
  if(ncol(binning_input) > 1){
    Dataset_3_bins[which(apply(binning_input_with_NA, 1, function(x) all(is.na(x)) )),-c(1)] <- NA
    #Dataset_3_bins[is.na(binning_input_with_NA)] <- NA
  }
  return(Dataset_3_bins)
}


data1 <- Three_bins("PXD010154_Ananth", "Tonsil")
data2 <- Three_bins("PXD010154_Ananth", "AdrenalGland")
data3 <- Three_bins("PXD010154_Ananth", "VermiformAppendix")
data4 <- Three_bins("PXD010154_Ananth", "BoneMarrow")
data5 <- Three_bins("PXD010154_Ananth", "Brain")
data6 <- Three_bins("PXD010154_Ananth", "Colon")
data7 <- Three_bins("PXD010154_Ananth", "Duodenum")
data8 <- Three_bins("PXD010154_Ananth", "UterineEndometrium")
data9 <- Three_bins("PXD010154_Ananth", "Esophagus")
data10 <- Three_bins("PXD010154_Ananth", "FallopianTubeOviduct")
data11 <- Three_bins("PXD010154_Ananth", "AdiposeTissue")
data12 <- Three_bins("PXD010154_Ananth", "GallBladder")
data13 <- Three_bins("PXD010154_Ananth", "Heart")
data14 <- Three_bins("PXD010154_Ananth", "Kidney")
data15 <- Three_bins("PXD010154_Ananth", "Liver")
data16 <- Three_bins("PXD010154_Ananth", "Lung")
data17 <- Three_bins("PXD010154_Ananth", "LymphNode")
data18 <- Three_bins("PXD010154_Ananth", "Ovary")
data19 <- Three_bins("PXD010154_Ananth", "Pancreas")
data20 <- Three_bins("PXD010154_Ananth", "PitutaryHypophysis")
data21 <- Three_bins("PXD010154_Ananth", "Placenta")
data22 <- Three_bins("PXD010154_Ananth", "Prostate")
data23 <- Three_bins("PXD010154_Ananth", "Rectum")
data24 <- Three_bins("PXD010154_Ananth", "SalivaryGland")
data25 <- Three_bins("PXD010154_Ananth", "SmallIntestine")
data26 <- Three_bins("PXD010154_Ananth", "SmoothMuscle")
data27 <- Three_bins("PXD010154_Ananth", "Spleen")
data28 <- Three_bins("PXD010154_Ananth", "Stomach")
data29 <- Three_bins("PXD010154_Ananth", "Testis")
data30 <- Three_bins("PXD010154_Ananth", "Thyroid")
data31 <- Three_bins("PXD010154_Ananth", "UrinaryBladder")
data32 <- Three_bins("PXD005819_33threads_yoda", "AnteriorPitutaryGland")
data33 <- Three_bins("PXD004143", "DorsoLateralPreFrontalCortex")
data34 <- Three_bins("PXD006233", "MiddleTemporalLobe")
data35 <- Three_bins("PXD012755", "CerebellarHemisphericCortex")
data36 <- Three_bins("PXD012755", "OccipitalCortex")
data37 <- Three_bins("PXD001608_30threads_yoda", "Colon")
data38 <- Three_bins("PXD002029", "Colon")
data39 <- Three_bins("PXD000547", "CorpusCallosum")
data40 <- Three_bins("PXD000548", "AnteriorTemporalLobe")
data41 <- Three_bins("PXD010271", "Pancreas")
data42 <- Three_bins("PXD010271", "Liver")
data43 <- Three_bins("PXD010271", "Ovary")
data44 <- Three_bins("PXD010271", "SubstantiaNigra")
data45 <- Three_bins("PXD004332", "PinealGland")
data46 <- Three_bins("PXD006675", "Aorta")
data47 <- Three_bins("PXD006675", "AorticValve")
data48 <- Three_bins("PXD006675", "LeftAtrium")
data49 <- Three_bins("PXD006675", "LeftVentricle")
data50 <- Three_bins("PXD006675", "MitralValve")
data51 <- Three_bins("PXD006675", "PulmonaryArtery")
data52 <- Three_bins("PXD006675", "PulmonaryValve")
data53 <- Three_bins("PXD006675", "PulmonaryVein")
data54 <- Three_bins("PXD006675", "RightAtrium")
data55 <- Three_bins("PXD006675", "RightVentricle")
data56 <- Three_bins("PXD006675", "AtrialSeptum")
data57 <- Three_bins("PXD006675", "VentricularSeptum")
data58 <- Three_bins("PXD006675", "TricuspidValve")
data59 <- Three_bins("PXD006675", "InferiorVenaCava")
data60 <- Three_bins("PXD008934", "LeftVentricle")
data61 <- Three_bins("Synapse-AD/ACT_DorsoLateralPreFrontalCortex", "DorsoLateralPreFrontalCortex")
data62 <- Three_bins("Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex", "DorsoLateralPreFrontalCortex")
data63 <- Three_bins("Synapse-AD/Banner_DorsoLateralPreFrontalCortex", "DorsoLateralPreFrontalCortex")
data64 <- Three_bins("Synapse-AD/BLSA_DorsoLateralPreFrontalCortex", "DorsoLateralPreFrontalCortex")
data65 <- Three_bins("Synapse-AD/BLSA_Precuneus", "Precuneus")
data66 <- Three_bins("Synapse-AD/Mayo_TemporalCortex", "TemporalCortex")
data67 <- Three_bins("Synapse-AD/MountSinai_FrontalPole", "FrontalPole")
data68 <- Three_bins("Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases", "DorsoLateralPreFrontalCortex")
data69 <- Three_bins("PXD012131", "Amygdala")
data70 <- Three_bins("PXD012131", "CaudateNucleus")
data71 <- Three_bins("PXD012131", "Cerebellum")
data72 <- Three_bins("PXD012131", "EntorhinalCortex")
data73 <- Three_bins("PXD012131", "InferiorParietalLobule")
data74 <- Three_bins("PXD012131", "FrontalGyrus")
data75 <- Three_bins("PXD012131", "NeoCortex")
data76 <- Three_bins("PXD012131", "SuperiorTemporalGyrus")
data77 <- Three_bins("PXD012131", "Thalamus")
data78 <- Three_bins("PXD012131", "VisualCortex")
data79 <- Three_bins("PXD020187", "UmblicalArtery")
data80 <- Three_bins("PXD015079", "PreFrontalCortex")
data81 <- Three_bins("PXD015079", "VermiformAppendix")


#data1 <- Three_bins("PXD010154_Ananth")
#data2 <- Three_bins("PXD005819_33threads_yoda")
#data3 <- Three_bins("PXD004143")
#data4 <- Three_bins("PXD006233")
#data5 <- Three_bins("PXD012755")
#data6 <- Three_bins("PXD001608_30threads_yoda")
#data7 <- Three_bins("PXD002029")
#data8 <- Three_bins("PXD000547")
#data9 <- Three_bins("PXD000548")
#data10 <- Three_bins("PXD010271")
#data11 <- Three_bins("PXD004332")
#data12 <- Three_bins("PXD006675")
#data13 <- Three_bins("PXD008934")
#data14 <- Three_bins("Synapse-AD/ACT_DorsoLateralPreFrontalCortex")
#data15 <- Three_bins("Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex")
#data16 <- Three_bins("Synapse-AD/Banner_DorsoLateralPreFrontalCortex")
#data17 <- Three_bins("Synapse-AD/BLSA_DorsoLateralPreFrontalCortex")
#data18 <- Three_bins("Synapse-AD/BLSA_Precuneus")
#data19 <- Three_bins("Synapse-AD/Mayo_TemporalCortex")
#data20 <- Three_bins("Synapse-AD/MountSinai_FrontalPole")
#data21 <- Three_bins("Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases")
#data22 <- Three_bins("PXD012131")
#data23 <- Three_bins("PXD020187")
#data24 <- Three_bins("PXD015079")



# Keep the minimum number of genes that are common in all datsets
Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene"), by.y=c("Gene"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(data1, data2)
merged_data <- Merge_data(merged_data, data3)
merged_data <- Merge_data(merged_data, data4)
merged_data <- Merge_data(merged_data, data5)
merged_data <- Merge_data(merged_data, data6)
merged_data <- Merge_data(merged_data, data7)
merged_data <- Merge_data(merged_data, data8)
merged_data <- Merge_data(merged_data, data9)
merged_data <- Merge_data(merged_data, data10)
merged_data <- Merge_data(merged_data, data11)
merged_data <- Merge_data(merged_data, data12)
merged_data <- Merge_data(merged_data, data13)
merged_data <- Merge_data(merged_data, data14)
merged_data <- Merge_data(merged_data, data15)
merged_data <- Merge_data(merged_data, data16)
merged_data <- Merge_data(merged_data, data17)
merged_data <- Merge_data(merged_data, data18)
merged_data <- Merge_data(merged_data, data19)
merged_data <- Merge_data(merged_data, data20)
merged_data <- Merge_data(merged_data, data21)
merged_data <- Merge_data(merged_data, data22)
merged_data <- Merge_data(merged_data, data23)
merged_data <- Merge_data(merged_data, data24)
merged_data <- Merge_data(merged_data, data25)
merged_data <- Merge_data(merged_data, data26)
merged_data <- Merge_data(merged_data, data27)
merged_data <- Merge_data(merged_data, data28)
merged_data <- Merge_data(merged_data, data29)
merged_data <- Merge_data(merged_data, data30)
merged_data <- Merge_data(merged_data, data31)
merged_data <- Merge_data(merged_data, data32)
merged_data <- Merge_data(merged_data, data33)
merged_data <- Merge_data(merged_data, data34)
merged_data <- Merge_data(merged_data, data35)
merged_data <- Merge_data(merged_data, data36)
merged_data <- Merge_data(merged_data, data37)
merged_data <- Merge_data(merged_data, data38)
merged_data <- Merge_data(merged_data, data39)
merged_data <- Merge_data(merged_data, data40)
merged_data <- Merge_data(merged_data, data41)
merged_data <- Merge_data(merged_data, data42)
merged_data <- Merge_data(merged_data, data43)
merged_data <- Merge_data(merged_data, data44)
merged_data <- Merge_data(merged_data, data45)
merged_data <- Merge_data(merged_data, data46)
merged_data <- Merge_data(merged_data, data47)
merged_data <- Merge_data(merged_data, data48)
merged_data <- Merge_data(merged_data, data49)
merged_data <- Merge_data(merged_data, data50)
merged_data <- Merge_data(merged_data, data51)
merged_data <- Merge_data(merged_data, data52)
merged_data <- Merge_data(merged_data, data53)
merged_data <- Merge_data(merged_data, data54)
merged_data <- Merge_data(merged_data, data55)
merged_data <- Merge_data(merged_data, data56)
merged_data <- Merge_data(merged_data, data57)
merged_data <- Merge_data(merged_data, data58)
merged_data <- Merge_data(merged_data, data59)
merged_data <- Merge_data(merged_data, data60)
merged_data <- Merge_data(merged_data, data61)
merged_data <- Merge_data(merged_data, data62)
merged_data <- Merge_data(merged_data, data63)
merged_data <- Merge_data(merged_data, data64)
merged_data <- Merge_data(merged_data, data65)
merged_data <- Merge_data(merged_data, data66)
merged_data <- Merge_data(merged_data, data67)
merged_data <- Merge_data(merged_data, data68)
merged_data <- Merge_data(merged_data, data69)
merged_data <- Merge_data(merged_data, data70)
merged_data <- Merge_data(merged_data, data71)
merged_data <- Merge_data(merged_data, data72)
merged_data <- Merge_data(merged_data, data73)
merged_data <- Merge_data(merged_data, data74)
merged_data <- Merge_data(merged_data, data75)
merged_data <- Merge_data(merged_data, data76)
merged_data <- Merge_data(merged_data, data77)
merged_data <- Merge_data(merged_data, data78)
merged_data <- Merge_data(merged_data, data79)
merged_data <- Merge_data(merged_data, data80)
merged_data <- Merge_data(merged_data, data81)

rm(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10,
   data11, data12, data13, data14, data15, data16, data17, data18, data19, data20,
   data21, data22, data23, data24, data25, data26, data27, data28, data29, data30,
   data31, data32, data33, data34, data35, data36, data37, data38, data39, data40,
   data41, data42, data43, data44, data45, data46, data47, data48, data49, data50,
   data51, data52, data53, data54, data55, data56, data57, data58, data59, data60,
   data61, data62, data63, data64, data65, data66, data67, data68, data69, data70,
   data71, data72, data73, data74, data75, data76, data77, data78, data79, data80, data81)

colnames(merged_data) <- gsub("PXD010154.Sample1.PitutaryHypophysis", "PXD010154.Sample1a.PitutaryHypophysis", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.Amygdala", "PXD012131.Sample1a.Amygdala", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.CaudateNucleus", "PXD012131.Sample1b.CaudateNucleus", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.Cerebellum", "PXD012131.Sample1c.Cerebellum", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.EntorhinalCortex", "PXD012131.Sample1e.EntorhinalCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.FrontalGyrus", "PXD012131.Sample1d.FrontalGyrus", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.InferiorParietalLobule", "PXD012131.Sample1f.InferiorParietalLobule", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample2.NeoCortex", "PXD012131.Sample2a.NeoCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample3.NeoCortex", "PXD012131.Sample3a.NeoCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample4.NeoCortex", "PXD012131.Sample4a.NeoCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.SuperiorTemporalGyrus", "PXD012131.Sample1g.SuperiorTemporalGyrus", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.Thalamus", "PXD012131.Sample1h.Thalamus", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012131.Sample1.VisualCortex", "PXD012131.Sample1i.VisualCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012755.Sample1.OccipitalCortex", "PXD012755.Sample1a.OccipitalCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012755.Sample2.OccipitalCortex", "PXD012755.Sample1b.OccipitalCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012755.Sample3.OccipitalCortex", "PXD012755.Sample1c.OccipitalCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012755.Sample4.OccipitalCortex", "PXD012755.Sample1d.OccipitalCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012755.Sample5.OccipitalCortex", "PXD012755.Sample1e.OccipitalCortex", colnames(merged_data))
colnames(merged_data) <- gsub("PXD012755.Sample6.OccipitalCortex", "PXD012755.Sample1f.OccipitalCortex", colnames(merged_data))

colnames(merged_data) <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                              "Brain", colnames(merged_data))

colnames(merged_data) <- gsub("PXD006675.Sample1.Aorta", "PXD006675.Sample1a.Aorta", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.Aorta", "PXD006675.Sample2a.Aorta", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.Aorta", "PXD006675.Sample3a.Aorta", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.AorticValve", "PXD006675.Sample1b.AorticValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.AorticValve", "PXD006675.Sample2b.AorticValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.AorticValve", "PXD006675.Sample3b.AorticValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.AtrialSeptum", "PXD006675.Sample1c.AtrialSeptum", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.AtrialSeptum", "PXD006675.Sample2c.AtrialSeptum", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.AtrialSeptum", "PXD006675.Sample3c.AtrialSeptum", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.InferiorVenaCava", "PXD006675.Sample1d.InferiorVenaCava", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.InferiorVenaCava", "PXD006675.Sample2d.InferiorVenaCava", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.InferiorVenaCava", "PXD006675.Sample3d.InferiorVenaCava", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.LeftAtrium", "PXD006675.Sample1e.LeftAtrium", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.LeftAtrium", "PXD006675.Sample2e.LeftAtrium", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.LeftAtrium", "PXD006675.Sample3e.LeftAtrium", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.LeftVentricle", "PXD006675.Sample1f.LeftVentricle", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.LeftVentricle", "PXD006675.Sample2f.LeftVentricle", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.LeftVentricle", "PXD006675.Sample3f.LeftVentricle", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.MitralValve", "PXD006675.Sample1g.MitralValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.MitralValve", "PXD006675.Sample2g.MitralValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.MitralValve", "PXD006675.Sample3g.MitralValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.PulmonaryArtery", "PXD006675.Sample1h.PulmonaryArtery", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.PulmonaryArtery", "PXD006675.Sample2h.PulmonaryArtery", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.PulmonaryArtery", "PXD006675.Sample3h.PulmonaryArtery", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.PulmonaryValve", "PXD006675.Sample1i.PulmonaryValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.PulmonaryValve", "PXD006675.Sample2i.PulmonaryValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.PulmonaryValve", "PXD006675.Sample3i.PulmonaryValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.PulmonaryVein", "PXD006675.Sample1j.PulmonaryVein", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.PulmonaryVein", "PXD006675.Sample2j.PulmonaryVein", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.PulmonaryVein", "PXD006675.Sample3j.PulmonaryVein", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.RightAtrium", "PXD006675.Sample1k.RightAtrium", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.RightAtrium", "PXD006675.Sample2k.RightAtrium", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.RightAtrium", "PXD006675.Sample3k.RightAtrium", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.RightVentricle", "PXD006675.Sample1l.RightVentricle", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.RightVentricle", "PXD006675.Sample2l.RightVentricle", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.RightVentricle", "PXD006675.Sample3l.RightVentricle", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.TricuspidValve", "PXD006675.Sample1m.TricuspidValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.TricuspidValve", "PXD006675.Sample2m.TricuspidValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.TricuspidValve", "PXD006675.Sample3m.TricuspidValve", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample1.VentricularSeptum", "PXD006675.Sample1n.VentricularSeptum", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample2.VentricularSeptum", "PXD006675.Sample2n.VentricularSeptum", colnames(merged_data))
colnames(merged_data) <- gsub("PXD006675.Sample3.VentricularSeptum", "PXD006675.Sample3n.VentricularSeptum", colnames(merged_data))


colnames(merged_data) <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                              "Heart", colnames(merged_data))

merged_data_long <- gather(merged_data, Datasets, Bins, colnames(merged_data[2]):colnames(merged_data[ncol(merged_data)]), factor_key=TRUE)
merged_data_long$Tissues <- gsub(".*\\.", "", merged_data_long$Datasets)

merged_data_long$Organs <- merged_data_long$Tissues
merged_data_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                                "Brain", merged_data_long$Organs)
merged_data_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                "Heart", merged_data_long$Organs)

#merged_data_long_aggregate <- aggregate(merged_data_long[,-c(2,4)], list("Gene" = merged_data_long$Gene, "Organs" = merged_data_long$Organs), median, na.rm =TRUE)
merged_data_long_aggregate <- aggregate(merged_data_long[,-c(2,4)], by=list(merged_data_long$Gene, merged_data_long$Organs), median, na.rm =TRUE)

merged_data_long_aggregate <- merged_data_long_aggregate[,-c(3,5)]
colnames(merged_data_long_aggregate) <- c("Gene", "Tissues", "Median_bins")

MQ_analysis_3bins<- spread(merged_data_long_aggregate, Tissues, Median_bins)

write.table(MQ_analysis_3bins, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/MQanalysis_3bins_all_samples.txt", sep = "\t", row.names = FALSE, quote = FALSE )


Gene_info <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Gene_distribution_in_organs-GeneNames-Median_intensities.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
Gene_info <- Gene_info[,c(1,2)]

MQ_analysis_3bins <- merge(x=Gene_info, y=MQ_analysis_3bins,
                           by.x=c("GeneName"), by.y=c("Gene"))


colnames(MQ_analysis_3bins) <- gsub("$", "_Manuscript_analysis", colnames(MQ_analysis_3bins), perl=TRUE)

MQ_Adipose <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$AdiposeTissue_Manuscript_analysis),c("GeneID_Manuscript_analysis", "AdiposeTissue_Manuscript_analysis"), drop=FALSE]
HPA_Adipose <- HPA_median_bins[!is.na(HPA_median_bins$adipose_tissue_HPA),c("Gene_HPA", "adipose_tissue_HPA"), drop=FALSE]
colnames(MQ_Adipose) <- c("Gene", "Level")
colnames(HPA_Adipose) <- c("Gene", "Level")
Total_Adipose <- unique(rbind(MQ_Adipose[,c("Gene"), drop=FALSE], HPA_Adipose[,c("Gene"), drop=FALSE]))
Common_Adipose      <- subset(MQ_Adipose, (Gene %in% HPA_Adipose$Gene))
Only_in_MQ_Adipose  <- subset(MQ_Adipose, !(Gene %in% HPA_Adipose$Gene))
Only_in_HPA_Adipose <- subset(HPA_Adipose, !(Gene %in% MQ_Adipose$Gene))

MQ_AdrenalGland <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$AdrenalGland_Manuscript_analysis),c("GeneID_Manuscript_analysis", "AdrenalGland_Manuscript_analysis"), drop=FALSE]
HPA_AdrenalGland <- HPA_median_bins[!is.na(HPA_median_bins$adrenal_gland_HPA),c("Gene_HPA","adrenal_gland_HPA"), drop=FALSE]
colnames(MQ_AdrenalGland) <- c("Gene", "Level")
colnames(HPA_AdrenalGland) <- c("Gene", "Level")
Total_AdrenalGland <- unique(rbind(MQ_AdrenalGland[,c("Gene"), drop=FALSE], HPA_AdrenalGland[,c("Gene"), drop=FALSE]))
Common_AdrenalGland      <- subset(MQ_AdrenalGland, (Gene %in% HPA_AdrenalGland$Gene))
Only_in_MQ_AdrenalGland  <- subset(MQ_AdrenalGland, !(Gene %in% HPA_AdrenalGland$Gene))
Only_in_HPA_AdrenalGland <- subset(HPA_AdrenalGland, !(Gene %in% MQ_AdrenalGland$Gene))

MQ_Appendix <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$VermiformAppendix_Manuscript_analysis),c("GeneID_Manuscript_analysis", "VermiformAppendix_Manuscript_analysis"), drop=FALSE]
HPA_Appendix <- HPA_median_bins[!is.na(HPA_median_bins$appendix_HPA),c("Gene_HPA","appendix_HPA"), drop=FALSE]
colnames(MQ_Appendix) <- c("Gene", "Level")
colnames(HPA_Appendix) <- c("Gene", "Level")
Total_Appendix <- unique(rbind(MQ_Appendix[,c("Gene"), drop=FALSE], HPA_Appendix[,c("Gene"), drop=FALSE]))
Common_Appendix      <- subset(MQ_Appendix, (Gene %in% HPA_Appendix$Gene))
Only_in_MQ_Appendix  <- subset(MQ_Appendix, !(Gene %in% HPA_Appendix$Gene))
Only_in_HPA_Appendix <- subset(HPA_Appendix, !(Gene %in% MQ_Appendix$Gene))

MQ_BoneMarrow <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$BoneMarrow_Manuscript_analysis),c("GeneID_Manuscript_analysis","BoneMarrow_Manuscript_analysis"), drop=FALSE]
HPA_BoneMarrow <- HPA_median_bins[!is.na(HPA_median_bins$bone_marrow_HPA),c("Gene_HPA","bone_marrow_HPA"), drop=FALSE]
colnames(MQ_BoneMarrow) <- c("Gene", "Level")
colnames(HPA_BoneMarrow) <- c("Gene", "Level")
Total_BoneMarrow <- unique(rbind(MQ_BoneMarrow[,c("Gene"), drop=FALSE], HPA_BoneMarrow[,c("Gene"), drop=FALSE]))
Common_BoneMarrow      <- subset(MQ_BoneMarrow, (Gene %in% HPA_BoneMarrow$Gene))
Only_in_MQ_BoneMarrow  <- subset(MQ_BoneMarrow, !(Gene %in% HPA_BoneMarrow$Gene))
Only_in_HPA_BoneMarrow <- subset(HPA_BoneMarrow, !(Gene %in% MQ_BoneMarrow$Gene))

MQ_Brain <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Brain_Manuscript_analysis),c("GeneID_Manuscript_analysis","Brain_Manuscript_analysis"), drop=FALSE]
HPA_Brain <- HPA_median_bins[!is.na(HPA_median_bins$Brain_HPA),c("Gene_HPA","Brain_HPA"), drop=FALSE]
colnames(MQ_Brain) <- c("Gene", "Level")
colnames(HPA_Brain) <- c("Gene", "Level")
Total_Brain <- unique(rbind(MQ_Brain[,c("Gene"), drop=FALSE], HPA_Brain[,c("Gene"), drop=FALSE]))
Common_Brain      <- subset(MQ_Brain, (Gene %in% HPA_Brain$Gene))
Only_in_MQ_Brain  <- subset(MQ_Brain, !(Gene %in% HPA_Brain$Gene))
Only_in_HPA_Brain <- subset(HPA_Brain, !(Gene %in% MQ_Brain$Gene))

#MQ_Breast <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Breast_Manuscript_analysis),c("GeneID_Manuscript_analysis","Breast_Manuscript_analysis"), drop=FALSE]
#HPA_Breast <- HPA_median_bins[!is.na(HPA_median_bins$breast_HPA),c("Gene_HPA","breast_HPA"), drop=FALSE]
#colnames(MQ_Breast) <- c("Gene", "Level")
#colnames(HPA_Breast) <- c("Gene", "Level")
#Total_Breast <- unique(rbind(MQ_Breast[,c("Gene"), drop=FALSE], HPA_Breast[,c("Gene"), drop=FALSE]))
#Common_Breast      <- subset(MQ_Breast, (Gene %in% HPA_Breast$Gene))
#Only_in_MQ_Breast  <- subset(MQ_Breast, !(Gene %in% HPA_Breast$Gene))
#Only_in_HPA_Breast <- subset(HPA_Breast, !(Gene %in% MQ_Breast$Gene))

MQ_Colon <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Colon_Manuscript_analysis),c("GeneID_Manuscript_analysis","Colon_Manuscript_analysis"), drop=FALSE]
HPA_Colon <- HPA_median_bins[!is.na(HPA_median_bins$colon_HPA),c("Gene_HPA","colon_HPA"), drop=FALSE]
colnames(MQ_Colon) <- c("Gene", "Level")
colnames(HPA_Colon) <- c("Gene", "Level")
Total_Colon <- unique(rbind(MQ_Colon[,c("Gene"), drop=FALSE], HPA_Colon[,c("Gene"), drop=FALSE]))
Common_Colon      <- subset(MQ_Colon, (Gene %in% HPA_Colon$Gene))
Only_in_MQ_Colon  <- subset(MQ_Colon, !(Gene %in% HPA_Colon$Gene))
Only_in_HPA_Colon <- subset(HPA_Colon, !(Gene %in% MQ_Colon$Gene))

MQ_Duodenum <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Duodenum_Manuscript_analysis),c("GeneID_Manuscript_analysis","Duodenum_Manuscript_analysis"), drop=FALSE]
HPA_Duodenum <- HPA_median_bins[!is.na(HPA_median_bins$duodenum_HPA),c("Gene_HPA","duodenum_HPA"), drop=FALSE]
colnames(MQ_Duodenum) <- c("Gene", "Level")
colnames(HPA_Duodenum) <- c("Gene", "Level")
Total_Duodenum <- unique(rbind(MQ_Duodenum[,c("Gene"), drop=FALSE], HPA_Duodenum[,c("Gene"), drop=FALSE]))
Common_Duodenum      <- subset(MQ_Duodenum, (Gene %in% HPA_Duodenum$Gene))
Only_in_MQ_Duodenum  <- subset(MQ_Duodenum, !(Gene %in% HPA_Duodenum$Gene))
Only_in_HPA_Duodenum <- subset(HPA_Duodenum, !(Gene %in% MQ_Duodenum$Gene))

MQ_Esophagus <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Esophagus_Manuscript_analysis),c("GeneID_Manuscript_analysis","Esophagus_Manuscript_analysis"), drop=FALSE]
HPA_Esophagus <- HPA_median_bins[!is.na(HPA_median_bins$esophagus_HPA),c("Gene_HPA","esophagus_HPA"), drop=FALSE]
colnames(MQ_Esophagus) <- c("Gene", "Level")
colnames(HPA_Esophagus) <- c("Gene", "Level")
Total_Esophagus <- unique(rbind(MQ_Esophagus[,c("Gene"), drop=FALSE], HPA_Esophagus[,c("Gene"), drop=FALSE]))
Common_Esophagus      <- subset(MQ_Esophagus, (Gene %in% HPA_Esophagus$Gene))
Only_in_MQ_Esophagus  <- subset(MQ_Esophagus, !(Gene %in% HPA_Esophagus$Gene))
Only_in_HPA_Esophagus <- subset(HPA_Esophagus, !(Gene %in% MQ_Esophagus$Gene))

MQ_FallopianTube <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$FallopianTubeOviduct_Manuscript_analysis),c("GeneID_Manuscript_analysis","FallopianTubeOviduct_Manuscript_analysis"), drop=FALSE]
HPA_FallopianTube <- HPA_median_bins[!is.na(HPA_median_bins$fallopian_tube_HPA),c("Gene_HPA","fallopian_tube_HPA"), drop=FALSE]
colnames(MQ_FallopianTube) <- c("Gene", "Level")
colnames(HPA_FallopianTube) <- c("Gene", "Level")
Total_FallopianTube <- unique(rbind(MQ_FallopianTube[,c("Gene"), drop=FALSE], HPA_FallopianTube[,c("Gene"), drop=FALSE]))
Common_FallopianTube      <- subset(MQ_FallopianTube, (Gene %in% HPA_FallopianTube$Gene))
Only_in_MQ_FallopianTube  <- subset(MQ_FallopianTube, !(Gene %in% HPA_FallopianTube$Gene))
Only_in_HPA_FallopianTube <- subset(HPA_FallopianTube, !(Gene %in% MQ_FallopianTube$Gene))

MQ_GallBladder <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$GallBladder_Manuscript_analysis),c("GeneID_Manuscript_analysis","GallBladder_Manuscript_analysis"), drop=FALSE]
HPA_GallBladder <- HPA_median_bins[!is.na(HPA_median_bins$gallbladder_HPA),c("Gene_HPA","gallbladder_HPA"), drop=FALSE]
colnames(MQ_GallBladder) <- c("Gene", "Level")
colnames(HPA_GallBladder) <- c("Gene", "Level")
Total_GallBladder <- unique(rbind(MQ_GallBladder[,c("Gene"), drop=FALSE], HPA_GallBladder[,c("Gene"), drop=FALSE]))
Common_GallBladder      <- subset(MQ_GallBladder, (Gene %in% HPA_GallBladder$Gene))
Only_in_MQ_GallBladder  <- subset(MQ_GallBladder, !(Gene %in% HPA_GallBladder$Gene))
Only_in_HPA_GallBladder <- subset(HPA_GallBladder, !(Gene %in% MQ_GallBladder$Gene))

MQ_Heart <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Heart_Manuscript_analysis),c("GeneID_Manuscript_analysis","Heart_Manuscript_analysis"), drop=FALSE]
HPA_Heart <- HPA_median_bins[!is.na(HPA_median_bins$heart_muscle_HPA),c("Gene_HPA","heart_muscle_HPA"), drop=FALSE]
colnames(MQ_Heart) <- c("Gene", "Level")
colnames(HPA_Heart) <- c("Gene", "Level")
Total_Heart <- unique(rbind(MQ_Heart[,c("Gene"), drop=FALSE], HPA_Heart[,c("Gene"), drop=FALSE]))
Common_Heart      <- subset(MQ_Heart, (Gene %in% HPA_Heart$Gene))
Only_in_MQ_Heart  <- subset(MQ_Heart, !(Gene %in% HPA_Heart$Gene))
Only_in_HPA_Heart <- subset(HPA_Heart, !(Gene %in% MQ_Heart$Gene))

MQ_Kidney <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Kidney_Manuscript_analysis),c("GeneID_Manuscript_analysis","Kidney_Manuscript_analysis"), drop=FALSE]
HPA_Kidney <- HPA_median_bins[!is.na(HPA_median_bins$kidney_HPA),c("Gene_HPA","kidney_HPA"), drop=FALSE]
colnames(MQ_Kidney) <- c("Gene", "Level")
colnames(HPA_Kidney) <- c("Gene", "Level")
Total_Kidney <- unique(rbind(MQ_Kidney[,c("Gene"), drop=FALSE], HPA_Kidney[,c("Gene"), drop=FALSE]))
Common_Kidney      <- subset(MQ_Kidney, (Gene %in% HPA_Kidney$Gene))
Only_in_MQ_Kidney  <- subset(MQ_Kidney, !(Gene %in% HPA_Kidney$Gene))
Only_in_HPA_Kidney <- subset(HPA_Kidney, !(Gene %in% MQ_Kidney$Gene))

MQ_Liver <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Liver_Manuscript_analysis),c("GeneID_Manuscript_analysis","Liver_Manuscript_analysis"), drop=FALSE]
HPA_Liver <- HPA_median_bins[!is.na(HPA_median_bins$liver_HPA),c("Gene_HPA","liver_HPA"), drop=FALSE]
colnames(MQ_Liver) <- c("Gene", "Level")
colnames(HPA_Liver) <- c("Gene", "Level")
Total_Liver <- unique(rbind(MQ_Liver[,c("Gene"), drop=FALSE], HPA_Liver[,c("Gene"), drop=FALSE]))
Common_Liver      <- subset(MQ_Liver, (Gene %in% HPA_Liver$Gene))
Only_in_MQ_Liver  <- subset(MQ_Liver, !(Gene %in% HPA_Liver$Gene))
Only_in_HPA_Liver <- subset(HPA_Liver, !(Gene %in% MQ_Liver$Gene))

MQ_Lung <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Lung_Manuscript_analysis),c("GeneID_Manuscript_analysis","Lung_Manuscript_analysis"), drop=FALSE]
HPA_Lung <- HPA_median_bins[!is.na(HPA_median_bins$lung_HPA),c("Gene_HPA","lung_HPA"), drop=FALSE]
colnames(MQ_Lung) <- c("Gene", "Level")
colnames(HPA_Lung) <- c("Gene", "Level")
Total_Lung <- unique(rbind(MQ_Lung[,c("Gene"), drop=FALSE], HPA_Lung[,c("Gene"), drop=FALSE]))
Common_Lung      <- subset(MQ_Lung, (Gene %in% HPA_Lung$Gene))
Only_in_MQ_Lung  <- subset(MQ_Lung, !(Gene %in% HPA_Lung$Gene))
Only_in_HPA_Lung <- subset(HPA_Lung, !(Gene %in% MQ_Lung$Gene))

MQ_LymphNode <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$LymphNode_Manuscript_analysis),c("GeneID_Manuscript_analysis","LymphNode_Manuscript_analysis"), drop=FALSE]
HPA_LymphNode <- HPA_median_bins[!is.na(HPA_median_bins$lymph_node_HPA),c("Gene_HPA","lymph_node_HPA"), drop=FALSE]
colnames(MQ_LymphNode) <- c("Gene", "Level")
colnames(HPA_LymphNode) <- c("Gene", "Level")
Total_LymphNode <- unique(rbind(MQ_LymphNode[,c("Gene"), drop=FALSE], HPA_LymphNode[,c("Gene"), drop=FALSE]))
Common_LymphNode      <- subset(MQ_LymphNode, (Gene %in% HPA_LymphNode$Gene))
Only_in_MQ_LymphNode  <- subset(MQ_LymphNode, !(Gene %in% HPA_LymphNode$Gene))
Only_in_HPA_LymphNode <- subset(HPA_LymphNode, !(Gene %in% MQ_LymphNode$Gene))

MQ_Ovary <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Ovary_Manuscript_analysis),c("GeneID_Manuscript_analysis","Ovary_Manuscript_analysis"), drop=FALSE]
HPA_Ovary <- HPA_median_bins[!is.na(HPA_median_bins$ovary_HPA),c("Gene_HPA","ovary_HPA"), drop=FALSE]
colnames(MQ_Ovary) <- c("Gene", "Level")
colnames(HPA_Ovary) <- c("Gene", "Level")
Total_Ovary <- unique(rbind(MQ_Ovary[,c("Gene"), drop=FALSE], HPA_Ovary[,c("Gene"), drop=FALSE]))
Common_Ovary      <- subset(MQ_Ovary, (Gene %in% HPA_Ovary$Gene))
Only_in_MQ_Ovary  <- subset(MQ_Ovary, !(Gene %in% HPA_Ovary$Gene))
Only_in_HPA_Ovary <- subset(HPA_Ovary, !(Gene %in% MQ_Ovary$Gene))

MQ_Pancreas <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Pancreas_Manuscript_analysis),c("GeneID_Manuscript_analysis","Pancreas_Manuscript_analysis"), drop=FALSE]
HPA_Pancreas <- HPA_median_bins[!is.na(HPA_median_bins$pancreas_HPA),c("Gene_HPA","pancreas_HPA"), drop=FALSE]
colnames(MQ_Pancreas) <- c("Gene", "Level")
colnames(HPA_Pancreas) <- c("Gene", "Level")
Total_Pancreas <- unique(rbind(MQ_Pancreas[,c("Gene"), drop=FALSE], HPA_Pancreas[,c("Gene"), drop=FALSE]))
Common_Pancreas      <- subset(MQ_Pancreas, (Gene %in% HPA_Pancreas$Gene))
Only_in_MQ_Pancreas  <- subset(MQ_Pancreas, !(Gene %in% HPA_Pancreas$Gene))
Only_in_HPA_Pancreas <- subset(HPA_Pancreas, !(Gene %in% MQ_Pancreas$Gene))

MQ_Placenta <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Placenta_Manuscript_analysis),c("GeneID_Manuscript_analysis","Placenta_Manuscript_analysis"), drop=FALSE]
HPA_Placenta <- HPA_median_bins[!is.na(HPA_median_bins$placenta_HPA),c("Gene_HPA","placenta_HPA"), drop=FALSE]
colnames(MQ_Placenta) <- c("Gene", "Level")
colnames(HPA_Placenta) <- c("Gene", "Level")
Total_Placenta <- unique(rbind(MQ_Placenta[,c("Gene"), drop=FALSE], HPA_Placenta[,c("Gene"), drop=FALSE]))
Common_Placenta      <- subset(MQ_Placenta, (Gene %in% HPA_Placenta$Gene))
Only_in_MQ_Placenta  <- subset(MQ_Placenta, !(Gene %in% HPA_Placenta$Gene))
Only_in_HPA_Placenta <- subset(HPA_Placenta, !(Gene %in% MQ_Placenta$Gene))

MQ_Prostate <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Prostate_Manuscript_analysis),c("GeneID_Manuscript_analysis","Prostate_Manuscript_analysis"), drop=FALSE]
HPA_Prostate <- HPA_median_bins[!is.na(HPA_median_bins$prostate_HPA),c("Gene_HPA","prostate_HPA"), drop=FALSE]
colnames(MQ_Prostate) <- c("Gene", "Level")
colnames(HPA_Prostate) <- c("Gene", "Level")
Total_Prostate <- unique(rbind(MQ_Prostate[,c("Gene"), drop=FALSE], HPA_Prostate[,c("Gene"), drop=FALSE]))
Common_Prostate      <- subset(MQ_Prostate, (Gene %in% HPA_Prostate$Gene))
Only_in_MQ_Prostate  <- subset(MQ_Prostate, !(Gene %in% HPA_Prostate$Gene))
Only_in_HPA_Prostate <- subset(HPA_Prostate, !(Gene %in% MQ_Prostate$Gene))

MQ_Rectum <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Rectum_Manuscript_analysis),c("GeneID_Manuscript_analysis","Rectum_Manuscript_analysis"), drop=FALSE]
HPA_Rectum <- HPA_median_bins[!is.na(HPA_median_bins$rectum_HPA),c("Gene_HPA","rectum_HPA"), drop=FALSE]
colnames(MQ_Rectum) <- c("Gene", "Level")
colnames(HPA_Rectum) <- c("Gene", "Level")
Total_Rectum <- unique(rbind(MQ_Rectum[,c("Gene"), drop=FALSE], HPA_Rectum[,c("Gene"), drop=FALSE]))
Common_Rectum      <- subset(MQ_Rectum, (Gene %in% HPA_Rectum$Gene))
Only_in_MQ_Rectum  <- subset(MQ_Rectum, !(Gene %in% HPA_Rectum$Gene))
Only_in_HPA_Rectum <- subset(HPA_Rectum, !(Gene %in% MQ_Rectum$Gene))

MQ_SalivaryGland <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$SalivaryGland_Manuscript_analysis),c("GeneID_Manuscript_analysis","SalivaryGland_Manuscript_analysis"), drop=FALSE]
HPA_SalivaryGland <- HPA_median_bins[!is.na(HPA_median_bins$salivary_gland_HPA),c("Gene_HPA","salivary_gland_HPA"), drop=FALSE]
colnames(MQ_SalivaryGland) <- c("Gene", "Level")
colnames(HPA_SalivaryGland) <- c("Gene", "Level")
Total_SalivaryGland <- unique(rbind(MQ_SalivaryGland[,c("Gene"), drop=FALSE], HPA_SalivaryGland[,c("Gene"), drop=FALSE]))
Common_SalivaryGland      <- subset(MQ_SalivaryGland, (Gene %in% HPA_SalivaryGland$Gene))
Only_in_MQ_SalivaryGland  <- subset(MQ_SalivaryGland, !(Gene %in% HPA_SalivaryGland$Gene))
Only_in_HPA_SalivaryGland <- subset(HPA_SalivaryGland, !(Gene %in% MQ_SalivaryGland$Gene))

MQ_SmallIntestine <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$SmallIntestine_Manuscript_analysis),c("GeneID_Manuscript_analysis","SmallIntestine_Manuscript_analysis"), drop=FALSE]
HPA_SmallIntestine <- HPA_median_bins[!is.na(HPA_median_bins$small_intestine_HPA),c("Gene_HPA","small_intestine_HPA"), drop=FALSE]
colnames(MQ_SmallIntestine) <- c("Gene", "Level")
colnames(HPA_SmallIntestine) <- c("Gene", "Level")
Total_SmallIntestine <- unique(rbind(MQ_SmallIntestine[,c("Gene"), drop=FALSE], HPA_SmallIntestine[,c("Gene"), drop=FALSE]))
Common_SmallIntestine      <- subset(MQ_SmallIntestine, (Gene %in% HPA_SmallIntestine$Gene))
Only_in_MQ_SmallIntestine  <- subset(MQ_SmallIntestine, !(Gene %in% HPA_SmallIntestine$Gene))
Only_in_HPA_SmallIntestine <- subset(HPA_SmallIntestine, !(Gene %in% MQ_SmallIntestine$Gene))

MQ_SmoothMuscle <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$SmoothMuscle_Manuscript_analysis),c("GeneID_Manuscript_analysis","SmoothMuscle_Manuscript_analysis"), drop=FALSE]
HPA_SmoothMuscle <- HPA_median_bins[!is.na(HPA_median_bins$smooth_muscle_HPA),c("Gene_HPA","smooth_muscle_HPA"), drop=FALSE]
colnames(MQ_SmoothMuscle) <- c("Gene", "Level")
colnames(HPA_SmoothMuscle) <- c("Gene", "Level")
Total_SmoothMuscle <- unique(rbind(MQ_SmoothMuscle[,c("Gene"), drop=FALSE], HPA_SmoothMuscle[,c("Gene"), drop=FALSE]))
Common_SmoothMuscle      <- subset(MQ_SmoothMuscle, (Gene %in% HPA_SmoothMuscle$Gene))
Only_in_MQ_SmoothMuscle  <- subset(MQ_SmoothMuscle, !(Gene %in% HPA_SmoothMuscle$Gene))
Only_in_HPA_SmoothMuscle <- subset(HPA_SmoothMuscle, !(Gene %in% MQ_SmoothMuscle$Gene))

MQ_Spleen <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Spleen_Manuscript_analysis),c("GeneID_Manuscript_analysis","Spleen_Manuscript_analysis"), drop=FALSE]
HPA_Spleen <- HPA_median_bins[!is.na(HPA_median_bins$spleen_HPA),c("Gene_HPA","spleen_HPA"), drop=FALSE]
colnames(MQ_Spleen) <- c("Gene", "Level")
colnames(HPA_Spleen) <- c("Gene", "Level")
Total_Spleen <- unique(rbind(MQ_Spleen[,c("Gene"), drop=FALSE], HPA_Spleen[,c("Gene"), drop=FALSE]))
Common_Spleen      <- subset(MQ_Spleen, (Gene %in% HPA_Spleen$Gene))
Only_in_MQ_Spleen  <- subset(MQ_Spleen, !(Gene %in% HPA_Spleen$Gene))
Only_in_HPA_Spleen <- subset(HPA_Spleen, !(Gene %in% MQ_Spleen$Gene))

MQ_Stomach <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Stomach_Manuscript_analysis),c("GeneID_Manuscript_analysis","Stomach_Manuscript_analysis"), drop=FALSE]
HPA_Stomach <- HPA_median_bins[!is.na(HPA_median_bins$Stomach_HPA),c("Gene_HPA","Stomach_HPA"), drop=FALSE]
colnames(MQ_Stomach) <- c("Gene", "Level")
colnames(HPA_Stomach) <- c("Gene", "Level")
Total_Stomach <- unique(rbind(MQ_Stomach[,c("Gene"), drop=FALSE], HPA_Stomach[,c("Gene"), drop=FALSE]))
Common_Stomach      <- subset(MQ_Stomach, (Gene %in% HPA_Stomach$Gene))
Only_in_MQ_Stomach  <- subset(MQ_Stomach, !(Gene %in% HPA_Stomach$Gene))
Only_in_HPA_Stomach <- subset(HPA_Stomach, !(Gene %in% MQ_Stomach$Gene))

MQ_Testis <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Testis_Manuscript_analysis),c("GeneID_Manuscript_analysis","Testis_Manuscript_analysis"), drop=FALSE]
HPA_Testis <- HPA_median_bins[!is.na(HPA_median_bins$testis_HPA),c("Gene_HPA","testis_HPA"), drop=FALSE]
colnames(MQ_Testis) <- c("Gene", "Level")
colnames(HPA_Testis) <- c("Gene", "Level")
Total_Testis <- unique(rbind(MQ_Testis[,c("Gene"), drop=FALSE], HPA_Testis[,c("Gene"), drop=FALSE]))
Common_Testis      <- subset(MQ_Testis, (Gene %in% HPA_Testis$Gene))
Only_in_MQ_Testis  <- subset(MQ_Testis, !(Gene %in% HPA_Testis$Gene))
Only_in_HPA_Testis <- subset(HPA_Testis, !(Gene %in% MQ_Testis$Gene))

MQ_Thyroid <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Thyroid_Manuscript_analysis),c("GeneID_Manuscript_analysis","Thyroid_Manuscript_analysis"), drop=FALSE]
HPA_Thyroid <- HPA_median_bins[!is.na(HPA_median_bins$thyroid_gland_HPA),c("Gene_HPA","thyroid_gland_HPA"), drop=FALSE]
colnames(MQ_Thyroid) <- c("Gene", "Level")
colnames(HPA_Thyroid) <- c("Gene", "Level")
Total_Thyroid <- unique(rbind(MQ_Thyroid[,c("Gene"), drop=FALSE], HPA_Thyroid[,c("Gene"), drop=FALSE]))
Common_Thyroid      <- subset(MQ_Thyroid, (Gene %in% HPA_Thyroid$Gene))
Only_in_MQ_Thyroid  <- subset(MQ_Thyroid, !(Gene %in% HPA_Thyroid$Gene))
Only_in_HPA_Thyroid <- subset(HPA_Thyroid, !(Gene %in% MQ_Thyroid$Gene))

MQ_Tonsil <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$Tonsil_Manuscript_analysis),c("GeneID_Manuscript_analysis","Tonsil_Manuscript_analysis"), drop=FALSE]
HPA_Tonsil <- HPA_median_bins[!is.na(HPA_median_bins$tonsil_HPA),c("Gene_HPA","tonsil_HPA"), drop=FALSE]
colnames(MQ_Tonsil) <- c("Gene", "Level")
colnames(HPA_Tonsil) <- c("Gene", "Level")
Total_Tonsil <- unique(rbind(MQ_Tonsil[,c("Gene"), drop=FALSE], HPA_Tonsil[,c("Gene"), drop=FALSE]))
Common_Tonsil      <- subset(MQ_Tonsil, (Gene %in% HPA_Tonsil$Gene))
Only_in_MQ_Tonsil  <- subset(MQ_Tonsil, !(Gene %in% HPA_Tonsil$Gene))
Only_in_HPA_Tonsil <- subset(HPA_Tonsil, !(Gene %in% MQ_Tonsil$Gene))

MQ_UrinaryBladder <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$UrinaryBladder_Manuscript_analysis),c("GeneID_Manuscript_analysis","UrinaryBladder_Manuscript_analysis"), drop=FALSE]
HPA_UrinaryBladder <- HPA_median_bins[!is.na(HPA_median_bins$urinary_bladder_HPA),c("Gene_HPA","urinary_bladder_HPA"), drop=FALSE]
colnames(MQ_UrinaryBladder) <- c("Gene", "Level")
colnames(HPA_UrinaryBladder) <- c("Gene", "Level")
Total_UrinaryBladder <- unique(rbind(MQ_UrinaryBladder[,c("Gene"), drop=FALSE], HPA_UrinaryBladder[,c("Gene"), drop=FALSE]))
Common_UrinaryBladder      <- subset(MQ_UrinaryBladder, (Gene %in% HPA_UrinaryBladder$Gene))
Only_in_MQ_UrinaryBladder  <- subset(MQ_UrinaryBladder, !(Gene %in% HPA_UrinaryBladder$Gene))
Only_in_HPA_UrinaryBladder <- subset(HPA_UrinaryBladder, !(Gene %in% MQ_UrinaryBladder$Gene))

MQ_UterineEndometrium <- MQ_analysis_3bins[!is.na(MQ_analysis_3bins$UterineEndometrium_Manuscript_analysis),c("GeneID_Manuscript_analysis","UterineEndometrium_Manuscript_analysis"), drop=FALSE]
HPA_UterineEndometrium <- HPA_median_bins[!is.na(HPA_median_bins$Endometrium_HPA),c("Gene_HPA","Endometrium_HPA"), drop=FALSE]
colnames(MQ_UterineEndometrium) <- c("Gene", "Level")
colnames(HPA_UterineEndometrium) <- c("Gene", "Level")
Total_UterineEndometrium <- unique(rbind(MQ_UterineEndometrium[,c("Gene"), drop=FALSE], HPA_UterineEndometrium[,c("Gene"), drop=FALSE]))
Common_UterineEndometrium      <- subset(MQ_UterineEndometrium, (Gene %in% HPA_UterineEndometrium$Gene))
Only_in_MQ_UterineEndometrium  <- subset(MQ_UterineEndometrium, !(Gene %in% HPA_UterineEndometrium$Gene))
Only_in_HPA_UterineEndometrium <- subset(HPA_UterineEndometrium, !(Gene %in% MQ_UterineEndometrium$Gene))


# table with percentage genes common between MQ analysis and HPA
Percent_common <- data.frame("Percent" = (nrow(Common_Adipose)/nrow(Total_Adipose))*100, "Organ" = "Adipose Tissue")
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_AdrenalGland)/nrow(Total_AdrenalGland))*100, "Organ" = "Adrenal Gland"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Appendix)/nrow(Total_Appendix))*100, "Organ" = "Vermiform Appendix"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_BoneMarrow)/nrow(Total_BoneMarrow))*100, "Organ" = "Bone Marrow"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Brain)/nrow(Total_Brain))*100, "Organ" = "Brain"))
#Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Breast)/nrow(Total_Breast))*100, "Organ" = "Breast"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Colon)/nrow(Total_Colon))*100, "Organ" = "Colon"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Duodenum)/nrow(Total_Duodenum))*100, "Organ" = "Duodenum"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Esophagus)/nrow(Total_Esophagus))*100, "Organ" = "Esophagus"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_FallopianTube)/nrow(Total_FallopianTube))*100, "Organ" = "Fallopian Tube Oviduct"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_GallBladder)/nrow(Total_GallBladder))*100, "Organ" = "Gallbladder"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Heart)/nrow(Total_Heart))*100, "Organ" = "Heart"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Kidney)/nrow(Total_Kidney))*100, "Organ" = "Kidney"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Liver)/nrow(Total_Liver))*100, "Organ" = "Liver"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Lung)/nrow(Total_Lung))*100, "Organ" = "Lung"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_LymphNode)/nrow(Total_LymphNode))*100, "Organ" = "Lymph Node"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Ovary)/nrow(Total_Ovary))*100, "Organ" = "Ovary"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Pancreas)/nrow(Total_Pancreas))*100, "Organ" = "Pancreas"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Placenta)/nrow(Total_Placenta))*100, "Organ" = "Placenta"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Prostate)/nrow(Total_Prostate))*100, "Organ" = "Prostate"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Rectum)/nrow(Total_Rectum))*100, "Organ" = "Rectum"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_SalivaryGland)/nrow(Total_SalivaryGland))*100, "Organ" = "Salivary Gland"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_SmallIntestine)/nrow(Total_SmallIntestine))*100, "Organ" = "Small Intestine"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_SmoothMuscle)/nrow(Total_SmoothMuscle))*100, "Organ" = "Smooth Muscle"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Spleen)/nrow(Total_Spleen))*100, "Organ" = "Spleen"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Stomach)/nrow(Total_Stomach))*100, "Organ" = "Stomach"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Testis)/nrow(Total_Testis))*100, "Organ" = "Testis"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Thyroid)/nrow(Total_Thyroid))*100, "Organ" = "Thyroid"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_Tonsil)/nrow(Total_Tonsil))*100, "Organ" = "Tonsil"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_UrinaryBladder)/nrow(Total_UrinaryBladder))*100, "Organ" = "Urinary Bladder"))
Percent_common <- rbind(Percent_common, data.frame("Percent" = (nrow(Common_UterineEndometrium)/nrow(Total_UterineEndometrium))*100, "Organ" = "Uterine Endometrium"))
Percent_common$Present <- rep("both in MQ analysis and HPA", nrow(Percent_common))

# table with percentage genes only identified in MQ analysis and NOT in HPA
Percent_only_in_MQ <- data.frame("Percent" = (nrow(Only_in_MQ_Adipose)/nrow(Total_Adipose))*100, "Organ" = "Adipose Tissue")
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_AdrenalGland)/nrow(Total_AdrenalGland))*100, "Organ" = "Adrenal Gland"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Appendix)/nrow(Total_Appendix))*100, "Organ" = "Vermiform Appendix"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_BoneMarrow)/nrow(Total_BoneMarrow))*100, "Organ" = "Bone Marrow"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Brain)/nrow(Total_Brain))*100, "Organ" = "Brain"))
#Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Breast)/nrow(Total_Breast))*100, "Organ" = "Breast"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Colon)/nrow(Total_Colon))*100, "Organ" = "Colon"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Duodenum)/nrow(Total_Duodenum))*100, "Organ" = "Duodenum"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Esophagus)/nrow(Total_Esophagus))*100, "Organ" = "Esophagus"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_FallopianTube)/nrow(Total_FallopianTube))*100, "Organ" = "Fallopian Tube Oviduct"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_GallBladder)/nrow(Total_GallBladder))*100, "Organ" = "Gallbladder"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Heart)/nrow(Total_Heart))*100, "Organ" = "Heart"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Kidney)/nrow(Total_Kidney))*100, "Organ" = "Kidney"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Liver)/nrow(Total_Liver))*100, "Organ" = "Liver"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Lung)/nrow(Total_Lung))*100, "Organ" = "Lung"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_LymphNode)/nrow(Total_LymphNode))*100, "Organ" = "Lymph Node"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Ovary)/nrow(Total_Ovary))*100, "Organ" = "Ovary"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Pancreas)/nrow(Total_Pancreas))*100, "Organ" = "Pancreas"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Placenta)/nrow(Total_Placenta))*100, "Organ" = "Placenta"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Prostate)/nrow(Total_Prostate))*100, "Organ" = "Prostate"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Rectum)/nrow(Total_Rectum))*100, "Organ" = "Rectum"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_SalivaryGland)/nrow(Total_SalivaryGland))*100, "Organ" = "Salivary Gland"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_SmallIntestine)/nrow(Total_SmallIntestine))*100, "Organ" = "Small Intestine"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_SmoothMuscle)/nrow(Total_SmoothMuscle))*100, "Organ" = "Smooth Muscle"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Spleen)/nrow(Total_Spleen))*100, "Organ" = "Spleen"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Stomach)/nrow(Total_Stomach))*100, "Organ" = "Stomach"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Testis)/nrow(Total_Testis))*100, "Organ" = "Testis"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Thyroid)/nrow(Total_Thyroid))*100, "Organ" = "Thyroid"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_Tonsil)/nrow(Total_Tonsil))*100, "Organ" = "Tonsil"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_UrinaryBladder)/nrow(Total_UrinaryBladder))*100, "Organ" = "Urinary Bladder"))
Percent_only_in_MQ <- rbind(Percent_only_in_MQ, data.frame("Percent" = (nrow(Only_in_MQ_UterineEndometrium)/nrow(Total_UterineEndometrium))*100, "Organ" = "Uterine Endometrium"))
Percent_only_in_MQ$Present <- rep("only in MQ analysis", nrow(Percent_only_in_MQ))

# table with percentage genes only identified in HPA and NOT in MQ analysis
Percent_only_in_HPA <- data.frame("Percent" = (nrow(Only_in_HPA_Adipose)/nrow(Total_Adipose))*100, "Organ" = "Adipose Tissue")
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_AdrenalGland)/nrow(Total_AdrenalGland))*100, "Organ" = "Adrenal Gland"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Appendix)/nrow(Total_Appendix))*100, "Organ" = "Vermiform Appendix"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_BoneMarrow)/nrow(Total_BoneMarrow))*100, "Organ" = "Bone Marrow"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Brain)/nrow(Total_Brain))*100, "Organ" = "Brain"))
#Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Breast)/nrow(Total_Breast))*100, "Organ" = "Breast"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Colon)/nrow(Total_Colon))*100, "Organ" = "Colon"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Duodenum)/nrow(Total_Duodenum))*100, "Organ" = "Duodenum"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Esophagus)/nrow(Total_Esophagus))*100, "Organ" = "Esophagus"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_FallopianTube)/nrow(Total_FallopianTube))*100, "Organ" = "Fallopian Tube Oviduct"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_GallBladder)/nrow(Total_GallBladder))*100, "Organ" = "Gallbladder"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Heart)/nrow(Total_Heart))*100, "Organ" = "Heart"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Kidney)/nrow(Total_Kidney))*100, "Organ" = "Kidney"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Liver)/nrow(Total_Liver))*100, "Organ" = "Liver"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Lung)/nrow(Total_Lung))*100, "Organ" = "Lung"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_LymphNode)/nrow(Total_LymphNode))*100, "Organ" = "Lymph Node"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Ovary)/nrow(Total_Ovary))*100, "Organ" = "Ovary"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Pancreas)/nrow(Total_Pancreas))*100, "Organ" = "Pancreas"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Placenta)/nrow(Total_Placenta))*100, "Organ" = "Placenta"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Prostate)/nrow(Total_Prostate))*100, "Organ" = "Prostate"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Rectum)/nrow(Total_Rectum))*100, "Organ" = "Rectum"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_SalivaryGland)/nrow(Total_SalivaryGland))*100, "Organ" = "Salivary Gland"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_SmallIntestine)/nrow(Total_SmallIntestine))*100, "Organ" = "Small Intestine"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_SmoothMuscle)/nrow(Total_SmoothMuscle))*100, "Organ" = "Smooth Muscle"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Spleen)/nrow(Total_Spleen))*100, "Organ" = "Spleen"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Stomach)/nrow(Total_Stomach))*100, "Organ" = "Stomach"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Testis)/nrow(Total_Testis))*100, "Organ" = "Testis"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Thyroid)/nrow(Total_Thyroid))*100, "Organ" = "Thyroid"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_Tonsil)/nrow(Total_Tonsil))*100, "Organ" = "Tonsil"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_UrinaryBladder)/nrow(Total_UrinaryBladder))*100, "Organ" = "Urinary Bladder"))
Percent_only_in_HPA <- rbind(Percent_only_in_HPA, data.frame("Percent" = (nrow(Only_in_HPA_UterineEndometrium)/nrow(Total_UterineEndometrium))*100, "Organ" = "Uterine Endometrium"))
Percent_only_in_HPA$Present <- rep("only in HPA", nrow(Percent_only_in_HPA))

plotdata <- rbind(Percent_common, Percent_only_in_MQ, Percent_only_in_HPA)

write.table(plotdata, "Comparison_of_genes_present_in_MQanalysis_and_HPA.txt", sep = "\t", row.names = FALSE, quote = FALSE )

plotdata <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Comparison_of_genes_present_in_MQanalysis_and_HPA.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
plotdata$Present <- gsub("MQ analysis", "this study", plotdata$Present, perl=TRUE)
  
ggplot(plotdata, aes(x=Organ, y=Percent, fill=factor(Present,levels=c("only in HPA","only in this study","both in this study and HPA")))) + 
  geom_bar(stat="identity") + 
  theme_bw()+
  labs(y="Percentage")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=12,),
        axis.title=element_text(size=12))+
  #scale_fill_manual(values = c("#440154FF","#2E6E8EFF","#FDE725FF"))  + 
  guides(fill=guide_legend(title="Present"))+
  ggtitle("Comparison of MQ data with genes identified from Human Protein Atlas (HPA)")


# Check correlation of expression

# First, just looking at genes being present (as 1) or absent (as 0)
Binary_correlation <- function(Tissue_MQ, Tissue_HPA){
  Merged_data <- merge(x = MQ_analysis_3bins[,c("GeneID_Manuscript_analysis", Tissue_MQ)], 
                       y = HPA_median_bins[,c("Gene_HPA", Tissue_HPA)],
                       by.x=c("GeneID_Manuscript_analysis"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
  
  Merged_data <- Merged_data[!is.na(Merged_data[[Tissue_MQ]]) | !is.na(Merged_data[[Tissue_HPA]]),]
  
  # undetected genes are given value 0 and detected genes have value 1 to help calculate r
  Merged_data[is.na(Merged_data)]  <- as.numeric(0)
  Merged_data[!is.na(Merged_data)] <- as.numeric(1)
  
  corel <- cor(Merged_data[[Tissue_MQ]], Merged_data[[Tissue_HPA]],  method = "pearson", use = "complete.obs")
  return(corel)
}

Present_absent_correl <- data.frame(Pearsons_r_for_TypeI=Binary_correlation("AdiposeTissue_Manuscript_analysis", "adipose_tissue_HPA"), Organ="Adipose Tissue")
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("AdrenalGland_Manuscript_analysis", "adrenal_gland_HPA"), Organ="Adrenal Gland"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("BoneMarrow_Manuscript_analysis", "bone_marrow_HPA"), Organ="Bone Marrow"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Brain_Manuscript_analysis", "Brain_HPA"), Organ="Brain"))
#Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Breast_Manuscript_analysis", "breast_HPA"), Organ="Breast"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Colon_Manuscript_analysis", "colon_HPA"), Organ="Colon"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Duodenum_Manuscript_analysis", "duodenum_HPA"), Organ="Duodenum"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Esophagus_Manuscript_analysis", "esophagus_HPA"), Organ="Esophagus"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("FallopianTubeOviduct_Manuscript_analysis", "fallopian_tube_HPA"), Organ="Fallopian Tube Oviduct"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("GallBladder_Manuscript_analysis", "gallbladder_HPA"), Organ="Gall Bladder"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Heart_Manuscript_analysis", "heart_muscle_HPA"), Organ="Heart"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Kidney_Manuscript_analysis", "kidney_HPA"), Organ="Kidney"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Liver_Manuscript_analysis", "liver_HPA"), Organ="Liver"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Lung_Manuscript_analysis", "lung_HPA"), Organ="Lung"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("LymphNode_Manuscript_analysis", "lymph_node_HPA"), Organ="Lymph Node"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Ovary_Manuscript_analysis", "ovary_HPA"), Organ="Ovary"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Pancreas_Manuscript_analysis", "pancreas_HPA"), Organ="Pancreas"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Placenta_Manuscript_analysis", "placenta_HPA"), Organ="Placenta"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Prostate_Manuscript_analysis", "prostate_HPA"), Organ="Prostate"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Rectum_Manuscript_analysis", "rectum_HPA"), Organ="Rectum"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("SalivaryGland_Manuscript_analysis", "salivary_gland_HPA"), Organ="Salivary Gland"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("SmallIntestine_Manuscript_analysis", "small_intestine_HPA"), Organ="Small Intestine"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("SmoothMuscle_Manuscript_analysis", "smooth_muscle_HPA"), Organ="Smooth Muscle"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Spleen_Manuscript_analysis", "spleen_HPA"), Organ="Spleen"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Stomach_Manuscript_analysis", "Stomach_HPA"), Organ="Stomach"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Testis_Manuscript_analysis", "testis_HPA"), Organ="Testis"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Thyroid_Manuscript_analysis", "thyroid_gland_HPA"), Organ="Thyroid"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("Tonsil_Manuscript_analysis", "tonsil_HPA"), Organ="Tonsil"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("UrinaryBladder_Manuscript_analysis", "urinary_bladder_HPA"), Organ="Urinary Bladder"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("UterineEndometrium_Manuscript_analysis", "Endometrium_HPA"), Organ="Uterine Endometrium"))
Present_absent_correl <- rbind(Present_absent_correl, data.frame(Pearsons_r_for_TypeI=Binary_correlation("VermiformAppendix_Manuscript_analysis", "appendix_HPA"), Organ="Vermiform Appendix"))
Present_absent_correl$Type1 <- rep("Detected(1)_or_undetected(0)", nrow(Present_absent_correl))

# Check against bin numbers 

Expression_correlation <- function(Tissue_MQ, Tissue_HPA){
  Merged_data <- merge(x = MQ_analysis_3bins[,c("GeneID_Manuscript_analysis", Tissue_MQ)], 
                       y = HPA_median_bins[,c("Gene_HPA", Tissue_HPA)],
                       by.x=c("GeneID_Manuscript_analysis"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)

  Merged_data <- Merged_data[!is.na(Merged_data[[Tissue_MQ]]) | !is.na(Merged_data[[Tissue_HPA]]),]

  #Merged_data[is.na(Merged_data)]  <- as.numeric(0)
  #Merged_data[!is.na(Merged_data)] <- as.numeric(1)

  corel <- cor(Merged_data[[Tissue_MQ]], Merged_data[[Tissue_HPA]],  method = "pearson", use = "complete.obs")
  return(corel)
}

Exp_correl <- data.frame(Pearsons_r_for_TypeII=Expression_correlation("AdiposeTissue_Manuscript_analysis", "adipose_tissue_HPA"), Organ="Adipose Tissue")
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("AdrenalGland_Manuscript_analysis", "adrenal_gland_HPA"), Organ="Adrenal Gland"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("BoneMarrow_Manuscript_analysis", "bone_marrow_HPA"), Organ="Bone Marrow"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Brain_Manuscript_analysis", "Brain_HPA"), Organ="Brain"))
#Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Breast_Manuscript_analysis", "breast_HPA"), Organ="Breast"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Colon_Manuscript_analysis", "colon_HPA"), Organ="Colon"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Duodenum_Manuscript_analysis", "duodenum_HPA"), Organ="Duodenum"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Esophagus_Manuscript_analysis", "esophagus_HPA"), Organ="Esophagus"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("FallopianTubeOviduct_Manuscript_analysis", "fallopian_tube_HPA"), Organ="Fallopian Tube Oviduct"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("GallBladder_Manuscript_analysis", "gallbladder_HPA"), Organ="Gall Bladder"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Heart_Manuscript_analysis", "heart_muscle_HPA"), Organ="Heart"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Kidney_Manuscript_analysis", "kidney_HPA"), Organ="Kidney"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Liver_Manuscript_analysis", "liver_HPA"), Organ="Liver"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Lung_Manuscript_analysis", "lung_HPA"), Organ="Lung"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("LymphNode_Manuscript_analysis", "lymph_node_HPA"), Organ="Lymph Node"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Ovary_Manuscript_analysis", "ovary_HPA"), Organ="Ovary"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Pancreas_Manuscript_analysis", "pancreas_HPA"), Organ="Pancreas"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Placenta_Manuscript_analysis", "placenta_HPA"), Organ="Placenta"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Prostate_Manuscript_analysis", "prostate_HPA"), Organ="Prostate"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Rectum_Manuscript_analysis", "rectum_HPA"), Organ="Rectum"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("SalivaryGland_Manuscript_analysis", "salivary_gland_HPA"), Organ="Salivary Gland"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("SmallIntestine_Manuscript_analysis", "small_intestine_HPA"), Organ="Small Intestine"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("SmoothMuscle_Manuscript_analysis", "smooth_muscle_HPA"), Organ="Smooth Muscle"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Spleen_Manuscript_analysis", "spleen_HPA"), Organ="Spleen"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Stomach_Manuscript_analysis", "Stomach_HPA"), Organ="Stomach"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Testis_Manuscript_analysis", "testis_HPA"), Organ="Testis"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Thyroid_Manuscript_analysis", "thyroid_gland_HPA"), Organ="Thyroid"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("Tonsil_Manuscript_analysis", "tonsil_HPA"), Organ="Tonsil"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("UrinaryBladder_Manuscript_analysis", "urinary_bladder_HPA"), Organ="Urinary Bladder"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("UterineEndometrium_Manuscript_analysis", "Endometrium_HPA"), Organ="Uterine Endometrium"))
Exp_correl <- rbind(Exp_correl, data.frame(Pearsons_r_for_TypeII=Expression_correlation("VermiformAppendix_Manuscript_analysis", "appendix_HPA"), Organ="Vermiform Appendix"))
Exp_correl$Type2 <- rep("Binned_values_against_HPA's_High(3)-Medium(2)-Low(1)", nrow(Exp_correl))


# Friedman's test

# Juan's suggestion: keep the current binning, and say 5-4 high, 2-3 medium and 1 low or something like that
#MQ_analysis_5bins_scaled_to_three <- MQ_analysis
#MQ_analysis_5bins_scaled_to_three[MQ_analysis_5bins_scaled_to_three >= 4 & MQ_analysis_5bins_scaled_to_three <= 5] <- 5
#MQ_analysis_5bins_scaled_to_three[MQ_analysis_5bins_scaled_to_three >= 2 & MQ_analysis_5bins_scaled_to_three <= 3] <- 3

#HPA_median_bins_scaled_to_three <- HPA_median_bins
#HPA_median_bins_scaled_to_three[HPA_median_bins_scaled_to_three >= 4 & HPA_median_bins_scaled_to_three <= 5] <- 3
#HPA_median_bins_scaled_to_three[HPA_median_bins_scaled_to_three >= 2 & HPA_median_bins_scaled_to_three <= 3] <- 2

Friedmans_correlation <- function(Tissue_MQ, Tissue_HPA){
  
  #Merged_data <- merge(x = MQ_analysis_3bins[,c("GeneID_Manuscript_analysis", Tissue_MQ)], 
  #                     y = HPA_median_bins[,c("Gene_HPA", Tissue_HPA)],
  #                     by.x=c("GeneID_Manuscript_analysis"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
  
  Merged_data <- merge(x = MQ_analysis_3bins[,c("Gene.ID", Tissue_MQ)], 
                       y = HPA_median_bins[,c("Gene_HPA", Tissue_HPA)],
                       by.x=c("Gene.ID"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
  
  
  Merged_data <- Merged_data[!is.na(Merged_data[[Tissue_MQ]]) | !is.na(Merged_data[[Tissue_HPA]]),]
  
  
  corel <- friedman.test(as.matrix(data.frame(Merged_data[[Tissue_MQ]], Merged_data[[Tissue_HPA]]),  ncol=2))
  return(corel$p.value)
}


Friedman_correl <- data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("AdiposeTissue", "adipose_tissue_HPA"), Organ="Adipose Tissue")
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("AdrenalGland", "adrenal_gland_HPA"), Organ="Adrenal Gland"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("BoneMarrow", "bone_marrow_HPA"), Organ="Bone Marrow"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Brain", "Brain_HPA"), Organ="Brain"))
#Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Breast", "breast_HPA"), Organ="Breast"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Colon", "colon_HPA"), Organ="Colon"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Duodenum", "duodenum_HPA"), Organ="Duodenum"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Esophagus", "esophagus_HPA"), Organ="Esophagus"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("FallopianTubeOviduct", "fallopian_tube_HPA"), Organ="Fallopian Tube Oviduct"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("GallBladder", "gallbladder_HPA"), Organ="Gall Bladder"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Heart", "heart_muscle_HPA"), Organ="Heart"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Kidney", "kidney_HPA"), Organ="Kidney"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Liver", "liver_HPA"), Organ="Liver"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Lung", "lung_HPA"), Organ="Lung"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("LymphNode", "lymph_node_HPA"), Organ="Lymph Node"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Ovary", "ovary_HPA"), Organ="Ovary"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Pancreas", "pancreas_HPA"), Organ="Pancreas"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Placenta", "placenta_HPA"), Organ="Placenta"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Prostate", "prostate_HPA"), Organ="Prostate"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Rectum", "rectum_HPA"), Organ="Rectum"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("SalivaryGland", "salivary_gland_HPA"), Organ="Salivary Gland"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("SmallIntestine", "small_intestine_HPA"), Organ="Small Intestine"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("SmoothMuscle", "smooth_muscle_HPA"), Organ="Smooth Muscle"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Spleen", "spleen_HPA"), Organ="Spleen"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Stomach", "Stomach_HPA"), Organ="Stomach"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Testis", "testis_HPA"), Organ="Testis"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Thyroid", "thyroid_gland_HPA"), Organ="Thyroid"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("Tonsil", "tonsil_HPA"), Organ="Tonsil"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("UrinaryBladder", "urinary_bladder_HPA"), Organ="Urinary Bladder"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("UterineEndometrium", "Endometrium_HPA"), Organ="Uterine Endometrium"))
Friedman_correl <- rbind(Friedman_correl, data.frame(Friedmans_pvalue_for_TypeIII=Friedmans_correlation("VermiformAppendix", "appendix_HPA"), Organ="Vermiform Appendix"))
Friedman_correl$Type3 <- rep("MQAnalysis_3bins", nrow(Friedman_correl))


Edit_distance <- function(Tissue_MQ, Tissue_HPA){
  
  Merged_data <- merge(x = MQ_analysis_3bins[,c("GeneID_Manuscript_analysis", Tissue_MQ)], 
                       y = HPA_median_bins[,c("Gene_HPA", Tissue_HPA)],
                       by.x=c("GeneID_Manuscript_analysis"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
  
  Merged_data <- Merged_data[!is.na(Merged_data[[Tissue_MQ]]) | !is.na(Merged_data[[Tissue_HPA]]),]
  
  Merged_data$absdiff_editdistance <- abs(Merged_data[[Tissue_MQ]]-Merged_data[[Tissue_HPA]])
  average_editdist <- mean(Merged_data$absdiff_editdistance, na.rm=TRUE)
  return(average_editdist)
}

Edit_dist <- data.frame(Mean_EditDistance=Edit_distance("AdiposeTissue_Manuscript_analysis", "adipose_tissue_HPA"), Organ="Adipose Tissue")
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("AdrenalGland_Manuscript_analysis", "adrenal_gland_HPA"), Organ="Adrenal Gland"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("BoneMarrow_Manuscript_analysis", "bone_marrow_HPA"), Organ="Bone Marrow"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Brain_Manuscript_analysis", "Brain_HPA"), Organ="Brain"))
#Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Breast_Manuscript_analysis", "breast_HPA"), Organ="Breast"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Colon_Manuscript_analysis", "colon_HPA"), Organ="Colon"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Duodenum_Manuscript_analysis", "duodenum_HPA"), Organ="Duodenum"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Esophagus_Manuscript_analysis", "esophagus_HPA"), Organ="Esophagus"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("FallopianTubeOviduct_Manuscript_analysis", "fallopian_tube_HPA"), Organ="Fallopian Tube Oviduct"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("GallBladder_Manuscript_analysis", "gallbladder_HPA"), Organ="Gall Bladder"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Heart_Manuscript_analysis", "heart_muscle_HPA"), Organ="Heart"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Kidney_Manuscript_analysis", "kidney_HPA"), Organ="Kidney"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Liver_Manuscript_analysis", "liver_HPA"), Organ="Liver"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Lung_Manuscript_analysis", "lung_HPA"), Organ="Lung"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("LymphNode_Manuscript_analysis", "lymph_node_HPA"), Organ="Lymph Node"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Ovary_Manuscript_analysis", "ovary_HPA"), Organ="Ovary"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Pancreas_Manuscript_analysis", "pancreas_HPA"), Organ="Pancreas"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Placenta_Manuscript_analysis", "placenta_HPA"), Organ="Placenta"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Prostate_Manuscript_analysis", "prostate_HPA"), Organ="Prostate"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Rectum_Manuscript_analysis", "rectum_HPA"), Organ="Rectum"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("SalivaryGland_Manuscript_analysis", "salivary_gland_HPA"), Organ="Salivary Gland"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("SmallIntestine_Manuscript_analysis", "small_intestine_HPA"), Organ="Small Intestine"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("SmoothMuscle_Manuscript_analysis", "smooth_muscle_HPA"), Organ="Smooth Muscle"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Spleen_Manuscript_analysis", "spleen_HPA"), Organ="Spleen"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Stomach_Manuscript_analysis", "Stomach_HPA"), Organ="Stomach"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Testis_Manuscript_analysis", "testis_HPA"), Organ="Testis"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Thyroid_Manuscript_analysis", "thyroid_gland_HPA"), Organ="Thyroid"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("Tonsil_Manuscript_analysis", "tonsil_HPA"), Organ="Tonsil"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("UrinaryBladder_Manuscript_analysis", "urinary_bladder_HPA"), Organ="Urinary Bladder"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("UterineEndometrium_Manuscript_analysis", "Endometrium_HPA"), Organ="Uterine Endometrium"))
Edit_dist <- rbind(Edit_dist, data.frame(Mean_EditDistance=Edit_distance("VermiformAppendix_Manuscript_analysis", "appendix_HPA"), Organ="Vermiform Appendix"))
Edit_dist$Type4 <- rep("Mean(abs(MQ_3bins-HPA_3bins))", nrow(Edit_dist))

#####Edit distances_matrix

### 1. ####################################
#All MQ analysis tissues against HPA tissues

#https://stackoverflow.com/questions/25384128/dividing-a-column-into-n-equal-groups-by-value
MQ_analysis_3bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/MCP/DDA_Manuscript_Supplementary_Files/Aggregate_Excel_rebinned_bins3.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


edit_distances_all <- data.frame()

for(i in 4:ncol(MQ_analysis_3bins)){
  for(j in 3:ncol(HPA_median_bins)){
      
#    for(i in 5:5){
#      for(j in 4:4){
      tmp <- data.frame()

      tmp <- merge(x = MQ_analysis_3bins[,c("Gene.ID", colnames(MQ_analysis_3bins)[i])] , 
                   y = HPA_median_bins[,c("Gene_HPA", colnames(HPA_median_bins)[j])] ,
                   by.x=c("Gene.ID"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
      
      tmp <- tmp[!is.na(tmp[[colnames(tmp)[2]]]) & !is.na(tmp[[colnames(tmp)[3]]]),]
      
      print(head(tmp))
      
      tmp$absdiff_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[3]]])
      average_editdist <- mean(tmp$absdiff_editdistance, na.rm=TRUE)
      
      result <- data.frame(Tissue1=colnames(tmp)[2], Tissue2=colnames(tmp)[3], average_absolute_editdistance=average_editdist)
      
      edit_distances_all <- rbind(edit_distances_all, result)
    }
}

#write.table(tmp, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance_AdrenalGland_only.txt", sep = "\t", row.names = TRUE, quote = FALSE )

  
edit_distances_all <- edit_distances_all[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix", edit_distances_all$Tissue1, ignore.case = TRUE), ]
edit_distances_all <- edit_distances_all[grep("adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_all$Tissue2, ignore.case =  TRUE), ]
edit_distances_all <- edit_distances_all[edit_distances_all$Tissue2 != "parathyroid_gland_HPA", ]


edit_distances_all$Tissue1 <- gsub("$", "_Prakash_etal_2021", edit_distances_all$Tissue1, perl=TRUE)
edit_distances_all$Tissue2 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("bone_marrow", "BoneMarrow", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("colon", "Colon", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("duodenum", "Duodenum", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("esophagus", "Esophagus", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("gallbladder", "GallBladder", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("heart_muscle", "Heart", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("kidney", "Kidney", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("liver", "Liver", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("lung", "Lung", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("lymph_node", "LymphNode", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("ovary", "Ovary", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("pancreas", "Pancreas", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("placenta", "Placenta", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("prostate", "Prostate", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("rectum", "Rectum", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("salivary_gland", "SalivaryGland", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("small_intestine", "SmallIntestine", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("spleen", "Spleen", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("testis", "Testis", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("tonsil", "Tonsil", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("Endometrium", "UterineEndometrium", edit_distances_all$Tissue2)
edit_distances_all$Tissue2 <- gsub("appendix", "VermiformAppendix", edit_distances_all$Tissue2)

edit_distance_matrix <- spread(edit_distances_all, Tissue2, average_absolute_editdistance)
rownames(edit_distance_matrix) <- edit_distance_matrix[,1]

edit_distance_matrix <- data.matrix(edit_distance_matrix[,-c(1)])

write.table(edit_distance_matrix, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance_matrix_Human_HPA.txt", sep = "\t", row.names = TRUE, quote = FALSE )

library(corrplot)

corrplot(edit_distance_matrix, is.corr = FALSE, type="full", tl.col="Black", tl.cex = 0.7, col = COL2('RdBu', 10))

par(oma=c(6,0,0,8)); 
heatmap.2(edit_distance_matrix, scale = "none", col = bluered(100), key.xlab="Edit distance",
          trace = "none", density.info = "none", cexRow=0.5, cexCol=0.75, margins = c(3,3))


# 2. ####################################
#All HPA tissues against HPA tissues
edit_distances_combined_MaxQuant_HPA  <- data.frame()

for(i in 3:ncol(HPA_median_bins)){
  for(j in 3:ncol(HPA_median_bins)){
    
    #    for(i in 5:5){
    #      for(j in 4:4){
    tmp <- data.frame()
    
    tmp <- merge(x = HPA_median_bins[,c("Gene_HPA", colnames(HPA_median_bins)[i])] , 
                 y = HPA_median_bins[,c("Gene_HPA", colnames(HPA_median_bins)[j])] ,
                 by.x=c("Gene_HPA"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
    
    tmp <- tmp[!is.na(tmp[[colnames(tmp)[2]]]) & !is.na(tmp[[colnames(tmp)[3]]]),]
    
    print(head(tmp))
    
    tmp$absdiff_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[3]]])
    average_editdist <- mean(tmp$absdiff_editdistance, na.rm=TRUE)
    
    result <- data.frame(Tissue1=colnames(tmp)[2], Tissue2=colnames(tmp)[3], average_absolute_editdistance=average_editdist)
    
    edit_distances_combined_MaxQuant_HPA  <- rbind(edit_distances_combined_MaxQuant_HPA , result)
  }
}

edit_distances_combined_MaxQuant_HPA  <- edit_distances_combined_MaxQuant_HPA [grep("adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_combined_MaxQuant_HPA $Tissue1, ignore.case =  TRUE), ]
edit_distances_combined_MaxQuant_HPA  <- edit_distances_combined_MaxQuant_HPA [grep("adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_combined_MaxQuant_HPA $Tissue2, ignore.case =  TRUE), ]
edit_distances_combined_MaxQuant_HPA  <- edit_distances_combined_MaxQuant_HPA [edit_distances_combined_MaxQuant_HPA $Tissue2 != "paraThyroidGland_HPA", ]
edit_distances_combined_MaxQuant_HPA  <- edit_distances_combined_MaxQuant_HPA [edit_distances_combined_MaxQuant_HPA $Tissue1 != "paraThyroidGland_HPA", ]


edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("\\.x", "", edit_distances_combined_MaxQuant_HPA $Tissue1, perl=TRUE)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("\\.y", "", edit_distances_combined_MaxQuant_HPA $Tissue2, perl=TRUE)

edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("bone_marrow", "BoneMarrow", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("colon", "Colon", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("duodenum", "Duodenum", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("esophagus", "Esophagus", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("gallbladder", "GallBladder", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("heart_muscle", "Heart", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("kidney", "Kidney", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("liver", "Liver", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("lung", "Lung", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("lymph_node", "LymphNode", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("ovary", "Ovary", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("pancreas", "Pancreas", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("placenta", "Placenta", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("prostate", "Prostate", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("rectum", "Rectum", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("salivary_gland", "SalivaryGland", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("small_intestine", "SmallIntestine", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("spleen", "Spleen", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("testis", "Testis", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("tonsil", "Tonsil", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("Endometrium", "UterineEndometrium", edit_distances_combined_MaxQuant_HPA $Tissue1)
edit_distances_combined_MaxQuant_HPA $Tissue1 <- gsub("appendix", "VermiformAppendix", edit_distances_combined_MaxQuant_HPA $Tissue1)

edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("bone_marrow", "BoneMarrow", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("colon", "Colon", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("duodenum", "Duodenum", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("esophagus", "Esophagus", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("gallbladder", "GallBladder", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("heart_muscle", "Heart", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("kidney", "Kidney", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("liver", "Liver", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("lung", "Lung", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("lymph_node", "LymphNode", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("ovary", "Ovary", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("pancreas", "Pancreas", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("placenta", "Placenta", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("prostate", "Prostate", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("rectum", "Rectum", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("salivary_gland", "SalivaryGland", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("small_intestine", "SmallIntestine", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("spleen", "Spleen", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("testis", "Testis", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("tonsil", "Tonsil", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("Endometrium", "UterineEndometrium", edit_distances_combined_MaxQuant_HPA $Tissue2)
edit_distances_combined_MaxQuant_HPA $Tissue2 <- gsub("appendix", "VermiformAppendix", edit_distances_combined_MaxQuant_HPA $Tissue2)

edit_distance_matrix_HPAs <- spread(edit_distances_combined_MaxQuant_HPA , Tissue2, average_absolute_editdistance)
rownames(edit_distance_matrix_HPAs) <- edit_distance_matrix_HPAs[,1]

edit_distance_matrix_HPAs <- data.matrix(edit_distance_matrix_HPAs[,-c(1)])

write.table(edit_distance_matrix_HPAs, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance_matrix_Human_HPA_vs_HPA.txt", sep = "\t", row.names = TRUE, quote = FALSE )

corrplot(edit_distance_matrix_HPAs, is.corr = FALSE, type="full", tl.col="Black", tl.cex = 0.7, col = COL2('RdBu', 10))

par(oma=c(6,0,0,8)); 
heatmap.2(edit_distance_matrix_HPAs, scale = "none", col = bluered(100), key.xlab="Edit distance",
          trace = "none", density.info = "none", cexRow=0.5, cexCol=0.75, margins = c(3,3))

# 3. ####################################
#All MaxQuant tissues against MaxQuant tissues

edit_distances_MaxQuants <- data.frame()

for(i in 4:ncol(MQ_analysis_3bins)){
  for(j in 4:ncol(MQ_analysis_3bins)){
    
    tmp <- data.frame()
    
    tmp <- merge(x = MQ_analysis_3bins[,c("Gene.ID", colnames(MQ_analysis_3bins)[i])] , 
                 y =  MQ_analysis_3bins[,c("Gene.ID", colnames(MQ_analysis_3bins)[j])] ,
                 by.x=c("Gene.ID"), by.y=c("Gene.ID"), all.x=FALSE, all.y=FALSE)
    
    tmp <- tmp[!is.na(tmp[[colnames(tmp)[2]]]) & !is.na(tmp[[colnames(tmp)[3]]]),]
    
    print(head(tmp))
    
    tmp$absdiff_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[3]]])
    average_editdist <- mean(tmp$absdiff_editdistance, na.rm=TRUE)
    
    result <- data.frame(Tissue1=colnames(tmp)[2], Tissue2=colnames(tmp)[3], average_absolute_editdistance=average_editdist)
    
    edit_distances_MaxQuants <- rbind(edit_distances_MaxQuants, result)
  }
}

edit_distances_MaxQuants <- edit_distances_MaxQuants[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix", edit_distances_MaxQuants$Tissue1, ignore.case = TRUE), ]
edit_distances_MaxQuants <- edit_distances_MaxQuants[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix", edit_distances_MaxQuants$Tissue2, ignore.case = TRUE), ]

edit_distances_MaxQuants$Tissue1 <- gsub("\\.x", "", edit_distances_MaxQuants$Tissue1, perl=TRUE)
edit_distances_MaxQuants$Tissue2 <- gsub("\\.y", "", edit_distances_MaxQuants$Tissue2, perl=TRUE)

edit_distances_MaxQuants$Tissue1 <- gsub("$", "_Prakash_etal_2021", edit_distances_MaxQuants$Tissue1, perl=TRUE)
edit_distances_MaxQuants$Tissue2 <- gsub("$", "_Prakash_etal_2021", edit_distances_MaxQuants$Tissue2, perl=TRUE)


edit_distance_matrix_MaxQuants <- spread(edit_distances_MaxQuants, Tissue2, average_absolute_editdistance)
rownames(edit_distance_matrix_MaxQuants) <- edit_distance_matrix_MaxQuants[,1]

edit_distance_matrix_MaxQuants <- data.matrix(edit_distance_matrix_MaxQuants[,-c(1)])

write.table(edit_distance_matrix_MaxQuants, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance_matrix_Human_MaxQuant_vs_MaxQuant.txt", sep = "\t", row.names = TRUE, quote = FALSE )

corrplot(edit_distance_matrix_MaxQuants, is.corr = FALSE, type="full", tl.col="Black", tl.cex = 0.7, col = COL2('RdBu', 10))

par(oma=c(10,0,0,8)); 
heatmap.2(edit_distance_matrix_MaxQuants, scale = "none", col = bluered(100), key.xlab="Edit distance",
          trace = "none", density.info = "none", cexRow=0.5, cexCol=0.75, margins = c(3,3))

#4. 
# MaxQuant+HPA vs MaxQuant+HPA
##############################


MQ_3bins_tmp <- MQ_analysis_3bins
colnames(MQ_3bins_tmp) <- gsub("$", "_Prakash_etal_2021", colnames(MQ_3bins_tmp), perl=TRUE)

combined_MaxQuant_HPA_bins <- merge(x = MQ_3bins_tmp, y = HPA_median_bins,
                                   by.x=c("Gene.ID_Prakash_etal_2021"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)

combined_MaxQuant_HPA_bins <- subset(combined_MaxQuant_HPA_bins, select=-c(Gene.name_HPA))

edit_distances_combined_MaxQuant_HPA <- data.frame()

for(i in 4:ncol(combined_MaxQuant_HPA_bins)){
  for(j in 4:ncol(combined_MaxQuant_HPA_bins)){
    
    tmp <- data.frame()
    
    tmp <- merge(x = combined_MaxQuant_HPA_bins[,c("Gene.ID_Prakash_etal_2021", colnames(combined_MaxQuant_HPA_bins)[i])] , 
                 y =  combined_MaxQuant_HPA_bins[,c("Gene.ID_Prakash_etal_2021", colnames(combined_MaxQuant_HPA_bins)[j])] ,
                 by.x=c("Gene.ID_Prakash_etal_2021"), by.y=c("Gene.ID_Prakash_etal_2021"), all.x=FALSE, all.y=FALSE)
    
    tmp <- tmp[!is.na(tmp[[colnames(tmp)[2]]]) & !is.na(tmp[[colnames(tmp)[3]]]),]
    
    print(head(tmp))
    
    tmp$absdiff_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[3]]])
    average_editdist <- mean(tmp$absdiff_editdistance, na.rm=TRUE)
    
    result <- data.frame(Tissue1=colnames(tmp)[2], Tissue2=colnames(tmp)[3], average_absolute_editdistance=average_editdist)
    
    edit_distances_combined_MaxQuant_HPA <- rbind(edit_distances_combined_MaxQuant_HPA, result)
  }
}

edit_distances_combined_MaxQuant_HPA <- edit_distances_combined_MaxQuant_HPA[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix|adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_combined_MaxQuant_HPA$Tissue1, ignore.case = TRUE), ]
edit_distances_combined_MaxQuant_HPA <- edit_distances_combined_MaxQuant_HPA[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix|adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_combined_MaxQuant_HPA$Tissue2, ignore.case =  TRUE), ]

edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("\\.x", "", edit_distances_combined_MaxQuant_HPA$Tissue1, perl=TRUE)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("\\.y", "", edit_distances_combined_MaxQuant_HPA$Tissue2, perl=TRUE)
edit_distances_combined_MaxQuant_HPA  <- edit_distances_combined_MaxQuant_HPA[edit_distances_combined_MaxQuant_HPA$Tissue1 != "parathyroid_gland_HPA", ]
edit_distances_combined_MaxQuant_HPA  <- edit_distances_combined_MaxQuant_HPA[edit_distances_combined_MaxQuant_HPA$Tissue2 != "parathyroid_gland_HPA", ]

edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("bone_marrow", "BoneMarrow", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("colon", "Colon", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("duodenum", "Duodenum", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("esophagus", "Esophagus", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("gallbladder", "GallBladder", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("heart_muscle", "Heart", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("kidney", "Kidney", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("liver", "Liver", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("lung", "Lung", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("lymph_node", "LymphNode", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("ovary", "Ovary", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("pancreas", "Pancreas", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("placenta", "Placenta", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("prostate", "Prostate", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("rectum", "Rectum", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("salivary_gland", "SalivaryGland", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("small_intestine", "SmallIntestine", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("spleen", "Spleen", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("testis", "Testis", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("tonsil", "Tonsil", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("Endometrium_HPA", "UterineEndometrium_HPA", edit_distances_combined_MaxQuant_HPA$Tissue1)
edit_distances_combined_MaxQuant_HPA$Tissue1 <- gsub("appendix", "VermiformAppendix", edit_distances_combined_MaxQuant_HPA$Tissue1)

edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("bone_marrow", "BoneMarrow", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("colon", "Colon", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("duodenum", "Duodenum", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("esophagus", "Esophagus", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("gallbladder", "GallBladder", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("heart_muscle", "Heart", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("kidney", "Kidney", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("liver", "Liver", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("lung", "Lung", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("lymph_node", "LymphNode", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("ovary", "Ovary", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("pancreas", "Pancreas", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("placenta", "Placenta", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("prostate", "Prostate", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("rectum", "Rectum", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("salivary_gland", "SalivaryGland", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("small_intestine", "SmallIntestine", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("spleen", "Spleen", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("testis", "Testis", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("tonsil", "Tonsil", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("Endometrium_HPA", "UterineEndometrium_HPA", edit_distances_combined_MaxQuant_HPA$Tissue2)
edit_distances_combined_MaxQuant_HPA$Tissue2 <- gsub("appendix", "VermiformAppendix", edit_distances_combined_MaxQuant_HPA$Tissue2)

edit_distances_combined_MaxQuant_HPA <- spread(edit_distances_combined_MaxQuant_HPA, Tissue2, average_absolute_editdistance)
rownames(edit_distances_combined_MaxQuant_HPA) <- edit_distances_combined_MaxQuant_HPA[,1]

edit_distances_combined_MaxQuant_HPA <- data.matrix(edit_distances_combined_MaxQuant_HPA[,-c(1)])

write.table(edit_distances_combined_MaxQuant_HPA, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance_matrix_Human_Alldatasets_vs_Alldatasets.txt", sep = "\t", row.names = TRUE, quote = FALSE )

corrplot(edit_distances_combined_MaxQuant_HPA, is.corr = FALSE, type="full", tl.col="Black", tl.cex = 0.7, col = COL2('RdBu', 10))

par(oma=c(10,0,0,8)); 
heatmap.2(edit_distances_combined_MaxQuant_HPA, scale = "none", col = bluered(100), key.xlab="Edit distance",
          trace = "none", density.info = "none", cexRow=0.5, cexCol=0.75, margins = c(3,3))

#heatmap(edit_distances_combined_MaxQuant_HPA, scale = "none")

#5. 
# MaxQuant+HPA vs MaxQuant+HPA (randomising a pair)
##############################

edit_distances_randomised_all <- data.frame()
all_values <- data.frame()

MQ_3bins_tmp <- MQ_analysis_3bins
colnames(MQ_3bins_tmp) <- gsub("$", "_Prakash_etal_2021", colnames(MQ_3bins_tmp), perl=TRUE)

combined_MaxQuant_HPA_bins <- merge(x = MQ_3bins_tmp, y = HPA_median_bins,
                                    by.x=c("Gene.ID_Prakash_etal_2021"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)

combined_MaxQuant_HPA_bins <- subset(combined_MaxQuant_HPA_bins, select=-c(Gene.name_HPA))


for(i in 4:ncol(combined_MaxQuant_HPA_bins)){
  for(j in 4:ncol(combined_MaxQuant_HPA_bins)){
    
    average_editdist <- c()
    average_randomised_editdist <- c()
    
    tmp <- data.frame()
    
    tmp <- merge(x = combined_MaxQuant_HPA_bins[,c("Gene.ID_Prakash_etal_2021", colnames(combined_MaxQuant_HPA_bins)[i])] , 
                 y =  combined_MaxQuant_HPA_bins[,c("Gene.ID_Prakash_etal_2021", colnames(combined_MaxQuant_HPA_bins)[j])] ,
                 by.x=c("Gene.ID_Prakash_etal_2021"), by.y=c("Gene.ID_Prakash_etal_2021"), all.x=FALSE, all.y=FALSE)
    
    tmp <- tmp[!is.na(tmp[[colnames(tmp)[2]]]) & !is.na(tmp[[colnames(tmp)[3]]]),]
    
    # sample does not work for rows less than 2
    if(nrow(tmp) > 1 ) {
      
      for(k in 1:10){
      #randomise HPA ten times (Andy's suggestion)
      print(head(tmp))
      tmp$randomised_Bin2 <- sample(tmp[,c(3)])
   
      # Correct edit distance: difference between Prakash_etal and HPA
      tmp$absdiff_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[3]]])
      average_editdist[k] <- mean(tmp$absdiff_editdistance, na.rm=TRUE)
    
      # Randomised edit distance: difference between Prakash_etal and randomised HPA
      tmp$absdiff_randomised_editdistance <- abs(tmp[[colnames(tmp)[2]]]-tmp[[colnames(tmp)[4]]])
      average_randomised_editdist[k] <- mean(tmp$absdiff_randomised_editdistance, na.rm=TRUE)
      } #end of for loop (k)
      
      # Difference between correct edit distance and randomised edit distance
      #tmp$editdistance_minus_randomised_distance <- abs(tmp$absdiff_editdistance-tmp$absdiff_randomised_editdistance)
      #average_editdist_minus_randomised <- mean(tmp$editdistance_minus_randomised_distance, na.rm=TRUE)
   
      #Andy suggested to use (random - true) instead of abs(true - random)
      #average_editdist_minus_randomised <- average_randomised_editdist-average_editdist
      average_editdist_minus_randomised <- mean(average_randomised_editdist)-mean(average_editdist)
      
      result <- data.frame(Tissue1=colnames(tmp)[2], Tissue2=colnames(tmp)[3], average_absolute_editdistance=mean(average_editdist), average_absolute_randomised_editdist=mean(average_randomised_editdist), editdist_minus_randomised=average_editdist_minus_randomised)
    
      edit_distances_randomised_all <- rbind(edit_distances_randomised_all, result)
     
      colnames(tmp)[1] <- "GeneID"

      Bin1 <- gsub(".*_Prakash","Prakash", colnames(tmp)[2])
      Bin1 <- gsub(".*_HPA","HPA", Bin1)
      Bin2 <- gsub(".*_Prakash","Prakash", colnames(tmp)[3])
      Bin2 <- gsub(".*_HPA","HPA", Bin2)
          
      tmp$Bin1_v._Bin2 <- paste(Bin1,Bin2, sep=".v.")
      tmp$Bin1_v._Bin2 <- gsub("\\.x","", tmp$Bin1_v._Bin2)
      tmp$Bin1_v._Bin2 <- gsub("\\.y","", tmp$Bin1_v._Bin2)

      Organ1 <- gsub("_.*","", colnames(tmp)[2])
      Organ2 <- gsub("_.*","", colnames(tmp)[3])
      tmp$Organ <- paste(Organ1, Organ2, sep=".v.")

      colnames(tmp)[2] <- gsub(".*_Prakash_etal_2021","Bin1", colnames(tmp)[2])
      colnames(tmp)[2] <- gsub(".*_HPA","Bin1", colnames(tmp)[2])
      colnames(tmp)[2] <- gsub("Bin1\\.x","Bin1", colnames(tmp)[2])
      colnames(tmp)[3] <- gsub(".*_Prakash_etal_2021","Bin2", colnames(tmp)[3])
      colnames(tmp)[3] <- gsub(".*_HPA","Bin2", colnames(tmp)[3])
      colnames(tmp)[3] <- gsub("Bin2\\.y","Bin2", colnames(tmp)[3])
      
      #print(head(tmp))
     # Saving all in a adataframe makes this loop run slow.
      #all_values <- rbind(all_values,tmp)
    }
  }
}

#all_values <- all_values[!grepl("breast|bronchus|cartilage|cervix|epididymis|Eye|hair|lactating|nasopharynx|oral|parathyroid|seminal|skeletal|soft|skin|thymus|vagina", all_values$Organ), ]
#write.table(all_values, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance__Human_randomised_AllValues_examples_new1.txt", sep = "\t", row.names = FALSE, quote = FALSE )

edit_distances_randomised_all <- edit_distances_randomised_all[,c("Tissue1", "Tissue2", "editdist_minus_randomised")]

edit_distances_randomised_all <- edit_distances_randomised_all[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix|adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_randomised_all$Tissue1, ignore.case = TRUE), ]
edit_distances_randomised_all <- edit_distances_randomised_all[grep("AdiposeTissue|AdrenalGland|BoneMarrow|Brain|Colon|Duodenum|Esophagus|FallopianTubeOviduct|GallBladder|Heart|Kidney|Liver|Lung|LymphNode|Ovary|Pancreas|Placenta|Prostate|Rectum|SalivaryGland|SmallIntestine|SmoothMuscle|Spleen|Stomach|Testis|Thyroid|Tonsil|UrinaryBladder|UterineEndometrium|VermiformAppendix|adipose_tissue|adrenal_gland|bone_marrow|Brain|colon|duodenum|esophagus|fallopian_tube|gallbladder|heart_muscle|kidney|liver|lung|lymph_node|ovary|pancreas|placenta|prostate|rectum|salivary_gland|small_intestine|smooth_muscle|spleen|Stomach|testis|thyroid_gland|tonsil|urinary_bladder|Endometrium|appendix", edit_distances_randomised_all$Tissue2, ignore.case =  TRUE), ]

edit_distances_randomised_all$Tissue1 <- gsub("\\.x", "", edit_distances_randomised_all$Tissue1, perl=TRUE)
edit_distances_randomised_all$Tissue2 <- gsub("\\.y", "", edit_distances_randomised_all$Tissue2, perl=TRUE)
edit_distances_randomised_all  <- edit_distances_randomised_all[edit_distances_randomised_all$Tissue1 != "parathyroid_gland_HPA", ]
edit_distances_randomised_all  <- edit_distances_randomised_all[edit_distances_randomised_all$Tissue2 != "parathyroid_gland_HPA", ]

edit_distances_randomised_all$Tissue1 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("bone_marrow", "BoneMarrow", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("colon", "Colon", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("duodenum", "Duodenum", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("esophagus", "Esophagus", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("gallbladder", "GallBladder", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("heart_muscle", "Heart", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("kidney", "Kidney", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("liver", "Liver", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("lung", "Lung", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("lymph_node", "LymphNode", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("ovary", "Ovary", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("pancreas", "Pancreas", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("placenta", "Placenta", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("prostate", "Prostate", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("rectum", "Rectum", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("salivary_gland", "SalivaryGland", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("small_intestine", "SmallIntestine", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("spleen", "Spleen", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("testis", "Testis", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("tonsil", "Tonsil", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("Endometrium_HPA", "UterineEndometrium_HPA", edit_distances_randomised_all$Tissue1)
edit_distances_randomised_all$Tissue1 <- gsub("appendix", "VermiformAppendix", edit_distances_randomised_all$Tissue1)

edit_distances_randomised_all$Tissue2 <- gsub("adipose_tissue", "AdiposeTissue", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("adrenal_gland", "AdrenalGland", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("bone_marrow", "BoneMarrow", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("colon", "Colon", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("duodenum", "Duodenum", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("esophagus", "Esophagus", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("fallopian_tube", "FallopianTubeOviduct", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("gallbladder", "GallBladder", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("heart_muscle", "Heart", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("kidney", "Kidney", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("liver", "Liver", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("lung", "Lung", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("lymph_node", "LymphNode", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("ovary", "Ovary", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("pancreas", "Pancreas", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("placenta", "Placenta", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("prostate", "Prostate", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("rectum", "Rectum", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("salivary_gland", "SalivaryGland", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("small_intestine", "SmallIntestine", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("smooth_muscle", "SmoothMuscle", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("spleen", "Spleen", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("testis", "Testis", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("thyroid_gland", "ThyroidGland", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("tonsil", "Tonsil", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("urinary_bladder", "UrinaryBladder", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("Endometrium_HPA", "UterineEndometrium_HPA", edit_distances_randomised_all$Tissue2)
edit_distances_randomised_all$Tissue2 <- gsub("appendix", "VermiformAppendix", edit_distances_randomised_all$Tissue2)

edit_distances_randomised_all <- spread(edit_distances_randomised_all, Tissue2, editdist_minus_randomised)
rownames(edit_distances_randomised_all) <- edit_distances_randomised_all[,1]

edit_distances_randomised_all <- data.matrix(edit_distances_randomised_all[,-c(1)])

colnames(edit_distances_randomised_all) <- gsub("Prakash_etal_2021", "This study", colnames(edit_distances_randomised_all))
rownames(edit_distances_randomised_all) <- gsub("Prakash_etal_2021", "This study", rownames(edit_distances_randomised_all))

write.table(edit_distances_randomised_all, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/HumanProteinAtlas/Edit_distance_matrix_Human_randomised_Alldatasets_new.txt", sep = "\t", row.names = TRUE, quote = FALSE )

par(oma=c(10,0,0,7)); 


heatmap.2(edit_distances_randomised_all, scale = "none", col = bluered(100), key.xlab="random edit distance difference",
          trace = "none", symbreaks = FALSE, density.info = "none", cexRow=0.75, cexCol=0.75, margins = c(3,3))

heatmap.2(edit_distances_randomised_all, scale = "none", col = bluered(100), key.xlab="random edit distance difference",
          trace = "none", symbreaks = FALSE, symm=F,symkey=F,density.info = "none", cexRow=0.5, cexCol=0.75, margins = c(3,3))



##################
#Spearman_correlation
##################
Spearman_correlation <- function(Tissue_MQ, Tissue_HPA){
  Merged_data <- merge(x = MQ_analysis_3bins[,c("Gene.ID", Tissue_MQ)], 
                       y = HPA_median_bins[,c("Gene_HPA", Tissue_HPA)],
                       by.x=c("Gene.ID"), by.y=c("Gene_HPA"), all.x=FALSE, all.y=FALSE)
  
  
  Merged_data <- Merged_data[!is.na(Merged_data[[Tissue_MQ]]) | !is.na(Merged_data[[Tissue_HPA]]),]
  
  
  res <- cor.test(x=Merged_data[[Tissue_MQ]], y=Merged_data[[Tissue_HPA]], method = 'spearman')
  return(res$estimate[[1]])
} 

Spearman_correl <- data.frame(Spearmans_rho=Spearman_correlation("AdiposeTissue", "adipose_tissue_HPA"), Organ="Adipose Tissue")
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("AdrenalGland", "adrenal_gland_HPA"), Organ="Adrenal Gland"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("BoneMarrow", "bone_marrow_HPA"), Organ="Bone Marrow"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Brain", "Brain_HPA"), Organ="Brain"))
#Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Breast", "breast_HPA"), Organ="Breast"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Colon", "colon_HPA"), Organ="Colon"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Duodenum", "duodenum_HPA"), Organ="Duodenum"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Esophagus", "esophagus_HPA"), Organ="Esophagus"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("FallopianTubeOviduct", "fallopian_tube_HPA"), Organ="Fallopian Tube Oviduct"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("GallBladder", "gallbladder_HPA"), Organ="Gall Bladder"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Heart", "heart_muscle_HPA"), Organ="Heart"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Kidney", "kidney_HPA"), Organ="Kidney"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Liver", "liver_HPA"), Organ="Liver"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Lung", "lung_HPA"), Organ="Lung"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("LymphNode", "lymph_node_HPA"), Organ="Lymph Node"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Ovary", "ovary_HPA"), Organ="Ovary"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Pancreas", "pancreas_HPA"), Organ="Pancreas"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Placenta", "placenta_HPA"), Organ="Placenta"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Prostate", "prostate_HPA"), Organ="Prostate"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Rectum", "rectum_HPA"), Organ="Rectum"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("SalivaryGland", "salivary_gland_HPA"), Organ="Salivary Gland"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("SmallIntestine", "small_intestine_HPA"), Organ="Small Intestine"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("SmoothMuscle", "smooth_muscle_HPA"), Organ="Smooth Muscle"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Spleen", "spleen_HPA"), Organ="Spleen"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Stomach", "Stomach_HPA"), Organ="Stomach"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Testis", "testis_HPA"), Organ="Testis"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Thyroid", "thyroid_gland_HPA"), Organ="Thyroid"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("Tonsil", "tonsil_HPA"), Organ="Tonsil"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("UrinaryBladder", "urinary_bladder_HPA"), Organ="Urinary Bladder"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("UterineEndometrium", "Endometrium_HPA"), Organ="Uterine Endometrium"))
Spearman_correl <- rbind(Spearman_correl, data.frame(Spearmans_rho=Spearman_correlation("VermiformAppendix", "appendix_HPA"), Organ="Vermiform Appendix"))
Spearman_correl$Type5 <- rep("MQAnalysis_3bins", nrow(Spearman_correl))



Correlations_HPA <- merge(x=Present_absent_correl, y=Exp_correl,
                          by.x=c("Organ"), by.y=c("Organ"), all.x=FALSE, all.y=FALSE)
Correlations_HPA <- merge(x=Correlations_HPA, y=Friedman_correl,
                          by.x=c("Organ"), by.y=c("Organ"), all.x=FALSE, all.y=FALSE)
Correlations_HPA <- merge(x=Correlations_HPA, y=Spearman_correl,
                          by.x=c("Organ"), by.y=c("Organ"), all.x=FALSE, all.y=FALSE)

Correlations_HPA <- merge(x=Correlations_HPA, y=Edit_dist,
                          by.x=c("Organ"), by.y=c("Organ"), all.x=FALSE, all.y=FALSE)

write.table(Correlations_HPA, "Correlation_with_HPA_data.txt", sep = "\t", row.names = FALSE, quote = FALSE )

#### Plot bin distribution of ids only seen in MQ and only seen in HPA

##### Proteins only in this study
All_organs_only_in_MQ <- rbind(Only_in_MQ_Adipose, Only_in_MQ_AdrenalGland, Only_in_MQ_Appendix, Only_in_MQ_BoneMarrow,
                               Only_in_MQ_Brain, Only_in_MQ_Colon, Only_in_MQ_Duodenum, Only_in_MQ_Esophagus, Only_in_MQ_FallopianTube,
                               Only_in_MQ_GallBladder, Only_in_MQ_Heart, Only_in_MQ_Kidney, Only_in_MQ_Liver,
                               Only_in_MQ_Lung, Only_in_MQ_LymphNode, Only_in_MQ_Ovary, Only_in_MQ_Pancreas,
                               Only_in_MQ_Placenta, Only_in_MQ_Prostate, Only_in_MQ_Rectum, Only_in_MQ_SalivaryGland,
                               Only_in_MQ_SmallIntestine, Only_in_MQ_SmoothMuscle, Only_in_MQ_Spleen, Only_in_MQ_Stomach,
                               Only_in_MQ_Testis, Only_in_MQ_Thyroid, Only_in_MQ_Tonsil, Only_in_MQ_UrinaryBladder,
                               Only_in_MQ_UterineEndometrium)
All_organs_only_in_MQ <- All_organs_only_in_MQ[,-c(2), drop=FALSE]
All_organs_only_in_MQ <- unique(All_organs_only_in_MQ)
All_organs_only_in_MQ <- merge(x=All_organs_only_in_MQ, y=MQ_analysis_3bins,
                               by.x=c("Gene"), by.y=c("Gene.ID"),
                               all.x=FALSE, all.y=FALSE)
All_organs_only_in_MQ_long <- gather(All_organs_only_in_MQ , Organs, Bins, colnames(All_organs_only_in_MQ[4]):colnames(All_organs_only_in_MQ[ncol(All_organs_only_in_MQ)]), factor_key=TRUE)
All_organs_only_in_MQ_long$Type <- rep("proteins identified only in this study", nrow(All_organs_only_in_MQ_long))

All_organs_only_in_MQ_percent <- All_organs_only_in_MQ_long[complete.cases(All_organs_only_in_MQ_long),] %>% 
                                  group_by(Bins, Type) %>% 
                                  summarise(total = n())

All_organs_only_in_MQ_percent <- group_by(All_organs_only_in_MQ_percent, Type) %>% 
                                  mutate(percent = (total/sum(total))*100)

All_organs_only_in_MQ_percent$Abundance_bins <- rep("from this study", nrow(All_organs_only_in_MQ_percent))

##### Proteins only in HPA
All_organs_only_in_HPA <- rbind(Only_in_HPA_Adipose, Only_in_HPA_AdrenalGland, Only_in_HPA_Appendix, Only_in_HPA_BoneMarrow,
                               Only_in_HPA_Brain, Only_in_HPA_Colon, Only_in_HPA_Duodenum, Only_in_HPA_Esophagus, Only_in_HPA_FallopianTube,
                               Only_in_HPA_GallBladder, Only_in_HPA_Heart, Only_in_HPA_Kidney, Only_in_HPA_Liver,
                               Only_in_HPA_Lung, Only_in_HPA_LymphNode, Only_in_HPA_Ovary, Only_in_HPA_Pancreas,
                               Only_in_HPA_Placenta, Only_in_HPA_Prostate, Only_in_HPA_Rectum, Only_in_HPA_SalivaryGland,
                               Only_in_HPA_SmallIntestine, Only_in_HPA_SmoothMuscle, Only_in_HPA_Spleen, Only_in_HPA_Stomach,
                               Only_in_HPA_Testis, Only_in_HPA_Thyroid, Only_in_HPA_Tonsil, Only_in_HPA_UrinaryBladder,
                               Only_in_HPA_UterineEndometrium)
All_organs_only_in_HPA <- All_organs_only_in_HPA[,-c(2), drop=FALSE]
All_organs_only_in_HPA <- unique(All_organs_only_in_HPA)
All_organs_only_in_HPA <- merge(x=All_organs_only_in_HPA, y=MQ_analysis_3bins,
                               by.x=c("Gene"), by.y=c("Gene.ID"),
                               all.x=FALSE, all.y=FALSE)
All_organs_only_in_HPA_long <- gather(All_organs_only_in_HPA , Organs, Bins, colnames(All_organs_only_in_HPA[4]):colnames(All_organs_only_in_HPA[ncol(All_organs_only_in_HPA)]), factor_key=TRUE)
All_organs_only_in_HPA_long$Type <- rep("proteins identified only in HPA", nrow(All_organs_only_in_HPA_long))

All_organs_only_in_HPA_percent <- All_organs_only_in_HPA_long[complete.cases(All_organs_only_in_HPA_long),] %>% 
  group_by(Bins, Type) %>% 
  summarise(total = n())

All_organs_only_in_HPA_percent <- group_by(All_organs_only_in_HPA_percent, Type) %>% 
  mutate(percent = (total/sum(total))*100)

All_organs_only_in_HPA_percent$Abundance_bins <- rep("from this study", nrow(All_organs_only_in_HPA_percent))


##### Proteins in both HPA and this study
All_organs_in_common <- rbind(Common_Adipose, Common_AdrenalGland, Common_Appendix, Common_BoneMarrow,
                              Common_Brain, Common_Colon, Common_Duodenum, Common_Esophagus, Common_FallopianTube,
                              Common_GallBladder, Common_Heart, Common_Kidney, Common_Liver,
                              Common_Lung, Common_LymphNode, Common_Ovary, Common_Pancreas,
                              Common_Placenta, Common_Prostate, Common_Rectum, Common_SalivaryGland,
                              Common_SmallIntestine, Common_SmoothMuscle, Common_Spleen, Common_Stomach,
                              Common_Testis, Common_Thyroid, Common_Tonsil, Common_UrinaryBladder,
                              Common_UterineEndometrium)
All_organs_in_common <- All_organs_in_common[,-c(2), drop=FALSE]
All_organs_in_common <- unique(All_organs_in_common)
All_organs_in_common <- merge(x=All_organs_in_common, y=MQ_analysis_3bins,
                              by.x=c("Gene"), by.y=c("Gene.ID"),
                              all.x=FALSE, all.y=FALSE)
All_organs_in_common_long <- gather(All_organs_in_common , Organs, Bins, colnames(All_organs_in_common[4]):colnames(All_organs_in_common[ncol(All_organs_in_common)]), factor_key=TRUE)
All_organs_in_common_long$Type <- rep("proteins identified in this study and HPA", nrow(All_organs_in_common_long))

All_organs_in_common_percent <- All_organs_in_common_long[complete.cases(All_organs_in_common_long),] %>% 
  group_by(Bins, Type) %>% 
  summarise(total = n())

All_organs_in_common_percent <- group_by(All_organs_in_common_percent, Type) %>% 
  mutate(percent = (total/sum(total))*100)

All_organs_in_common_percent$Abundance_bins <- rep("from this study", nrow(All_organs_in_common_percent))


#Getting bin data as in HPA for genes identified in HPA and in this study 

HPA_bins_of_common_proteins <- merge(x=tmp_subdata_aggregate, y=All_organs_in_common[,c("Gene", "Gene.Name")],
                                     by.x=c("Gene","Gene.name"), by.y=c("Gene","Gene.Name"),
                                     all.x=FALSE, all.y=FALSE) 

HPA_bins_of_common_proteins <- HPA_bins_of_common_proteins[complete.cases(HPA_bins_of_common_proteins),] %>% 
  group_by(Level) %>% 
  summarise(total = n())
HPA_bins_of_common_proteins$percent <- (HPA_bins_of_common_proteins$total/sum(HPA_bins_of_common_proteins$total))*100
HPA_bins_of_common_proteins$Type <- rep("proteins identified in this study and HPA", nrow(HPA_bins_of_common_proteins))
colnames(HPA_bins_of_common_proteins)[1] <- "Bins"

HPA_bins_of_common_proteins <- HPA_bins_of_common_proteins[,c("Bins","Type","total","percent")]

HPA_bins_of_common_proteins$Abundance_bins <- rep("from HPA", nrow(HPA_bins_of_common_proteins))


plotdata <- rbind(All_organs_only_in_MQ_percent, All_organs_only_in_HPA_percent, All_organs_in_common_percent, HPA_bins_of_common_proteins)

ggplot(plotdata, aes(fill=as.factor(Abundance_bins), x=Bins, y=percent)) + 
  geom_col(colour="grey",position = position_dodge(width = 0.25))+ 
  xlab("Protein abundance bins")+
  ylab("Percent")+
  theme_bw()+
  labs(fill="Abundance bins sourced")+
  #scale_fill_manual(values = c("#A3A500", "#00B0F6"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = c(0.85, 0.8))+
  ggtitle("Comparison of proteins identified between this study and Human Protein Atlas (HPA)")+
  facet_wrap(~Type)





#### Calculating bin value percentages as downloaded from HPA
all_HPA_percent <- tmp_subdata_aggregate[complete.cases(tmp_subdata_aggregate),] %>% 
  group_by(Level) %>% 
  summarise(total = n())

all_HPA_percent$percent <- (all_HPA_percent$total/sum(all_HPA_percent$total))*100
all_HPA_percent$Type <- rep("HPA", nrow(all_HPA_percent))
colnames(all_HPA_percent)[1] <- "Bins"

all_MaxQuant_percent <- gather(MQ_analysis_3bins, Organs, Bins, colnames(MQ_analysis_3bins[4]):colnames(MQ_analysis_3bins[ncol(MQ_analysis_3bins)]), factor_key=TRUE)
all_MaxQuant_percent <- all_MaxQuant_percent[complete.cases(all_MaxQuant_percent),] %>% 
  group_by(Bins) %>% 
  summarise(total = n())

all_MaxQuant_percent$percent <- (all_MaxQuant_percent$total/sum(all_MaxQuant_percent$total))*100
all_MaxQuant_percent$Type <- rep("This study", nrow(all_MaxQuant_percent))


all_percents_plotdata <- rbind(all_HPA_percent, all_MaxQuant_percent)

ggplot(all_percents_plotdata , aes(fill=as.factor(Type), x=Bins, y=percent)) + 
  geom_col(colour="grey", position = position_dodge(width = 0.25))+ 
  xlab("Protein abundance bins")+
  ylab("Percent")+
  theme_bw()+
  labs(fill="Bins")+
  #scale_fill_manual(values = c("#A3A500", "#00B0F6"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #theme(legend.position = "none")+
  ggtitle("Background distribution of all binned data in this study and\nHuman Protein Atlas (HPA)")
