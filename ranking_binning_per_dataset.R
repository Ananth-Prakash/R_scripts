library(ggplot2)
library(dplyr)
library(stats)
library(tidyr)
library(viridis)
library(grid)
library(RColorBrewer)
library(heatmap3)
library(stringr)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
library(sva)

BiocManager::install("limma")
library(limma)


# Input: ppb (FOT) normalised intensity values from MaxQuant,
#        Normal tissue samples only
#        Rename coloumns to add PXDID and Sample number
# 

# breast 2 : PXD001325, PXD012431
# heart  3 : PXD006675 (sub regions), PXD008934, PXD008722 (sub regions)
# brain    : PXD005819, PXD004143, PXD006233, PXD012755 (sub regions), PXD000547, PXD000548, PXD010271, PXD004332
#          : Synapse-ACT-DLPFC (normal), Synapse-Aging-DLPFC (normal), Synapse-Banner-DLPFC (normal) 
# colon    : PXD001608, PXD002029

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/')

PXDID <- "PXD001608"

tmp  <- read.table(file = paste("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/", PXDID, "/proteinGroups_ppb_final-tissue_names.txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

#tmp  <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases/proteinGroups_ppb_final-tissue_names.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")



#Some of the gene entries (ex. IGHA2 has two Ensembl gene ids ENSG00000211890 & ENSG00000276173
#                          ex. IGHV2-70 has two Ensembl gene ids ENSG00000274576 & ENSG00000282453)
# because of this there are duplicate entries of such genes. These are aggregated by taking the median of them

colnames(tmp)[2] <- "Gene.Symbol"
dataset <- aggregate(tmp[,-c(1,2)], list(tmp$Gene.Symbol), median, na.rm =TRUE)
#dataset <- subset(dataset, select = -c(Gene.Symbol))
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


Binning <- function(data_input, gene_names) {
  # Make a data frame of sample + protein/peptide abundances. 
  data = data.frame(data_input)
  
  # Rank data.
  data_ranked = as.data.frame(t(convert_to_ranked_data(data)))
  colnames(data_ranked) <- gene_names
  
  # Different size of bins to try.
  binning_size_approaches = list(
    bins_2 = 1,
    bins_5 = 4,
    bins_10 = 9,
    bins_100 = 99,
    bins_1000 = 999,
    bins_10000 = 9999
  )
  
  # Transform data using different bin sizes.
  results = lapply(binning_size_approaches, function(binning_size_approach) {
    # 'Bin' data by scaling the ranked data.
    binned_data = as.data.frame(t(convert_to_scaled_data(data_ranked, binning_size_approach)))
  })
}

# If a dataset has different tissues or tissue regions, do not bin all together,
# these tissues/regions must be separated and binned separately.

colnames(dataset)

binning_input <- dataset[,-c(1)]
#binning_input <- dataset[,c(2:9), drop=FALSE]

#keep a copy of the input before changing NAs to 0
binning_input_with_NA <- binning_input

# non-detected genes (NA or 0 values) are assigned 0 values here to help in binning.
binning_input[is.na(binning_input)] <- 0


Dataset_all_bins <- Binning(binning_input, gene_names)


# Consider only data grouped into 5 bins
Get_Bin_5 <- function(binned_data_list){
  
  Bins_5 <- as.data.frame(t(binned_data_list$bins_5))
  # rename bins 0 to 4 --> 1 to 5
  # lowest intensities are in bin1 and highest in bin5
  Bins_5[Bins_5 == 4] <- 5
  Bins_5[Bins_5 == 3] <- 4
  Bins_5[Bins_5 == 2] <- 3
  Bins_5[Bins_5 == 1] <- 2
  Bins_5[Bins_5 == 0] <- 1
  
  Bins_5 <- tibble::rownames_to_column(Bins_5, "Gene")
}

Dataset_5_bins <- Get_Bin_5(Dataset_all_bins)


# Now reassign NA values back to those genes which were undetected or 0 in ALL samples 
# NOTE: Only those genes are reassigned back to NA that have missing values in
# all samples (ie., if at least one sample had a value then it is left as 0)
# this particular reassigning is done after the binning (later in the code).


Dataset_5_bins[which(apply(binning_input_with_NA, 1, function(x) all(is.na(x)) )),-c(1)] <- NA

write.table(Dataset_5_bins, file = paste("Binned_expression_", PXDID ,"_Amygdala", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )
#write.table(Dataset_5_bins, file = paste("Binned_expression_syn21443008_Brain_UniPenn-DLPFC.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )


# Count the number of occurrences of each bin
Bins_count <- function(dataset){
  Bin_1 <- apply(dataset, 1, function(x) length(which(x == 1)) )
  Bin_2 <- apply(dataset, 1, function(x) length(which(x == 2)) )
  Bin_3 <- apply(dataset, 1, function(x) length(which(x == 3)) )
  Bin_4 <- apply(dataset, 1, function(x) length(which(x == 4)) )
  Bin_5 <- apply(dataset, 1, function(x) length(which(x == 5)) )
  
  result <- as.data.frame(cbind(Bin_1, Bin_2, Bin_3, Bin_4, Bin_5))
}

Bincounts <- cbind( Gene=Dataset_5_bins[,c(1)], as.data.frame( Bins_count ( as.data.frame(Dataset_5_bins[,-c(1)])  ) ) )
write.table(Bincounts, file = paste("Bin_counts_", PXDID, "_Amygdala", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )
#write.table(Bincounts, file = paste("Bin_counts_syn21443008_Brain_UniPenn-DLPFC.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

# After this when comparing different datasets of the same tissue type, merge them
# this way there will be same genes between datasets when comparing

##### Comparing datasets

##### BRAIN

dataset1 <- read.table("Binned_expression_PXD000547_Brain-CorpusCallosum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table("Binned_expression_PXD000548_Brain-AnteriorTemporalLobe.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset3 <- read.table("Binned_expression_PXD004143_Brain-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table("Binned_expression_PXD004332_Brain-PinealGland.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table("Binned_expression_PXD005819_Brain-AnteriorPituitaryGland.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table("Binned_expression_PXD006233_Brain-MiddleTemporalLobe.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table("Binned_expression_PXD010154_Brain.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table("Binned_expression_PXD010154_Brain-PitutaryHypophysis.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table("Binned_expression_PXD010271_Brain-SubstantiaNigra.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset10 <- read.table("Binned_expression_PXD012131_Brain-Amygdala.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset11 <- read.table("Binned_expression_PXD012131_Brain-CaudateNucleus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset12 <- read.table("Binned_expression_PXD012131_Brain-Cerebellum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset13 <- read.table("Binned_expression_PXD012131_Brain-EntorhinalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset14 <- read.table("Binned_expression_PXD012131_Brain-FrontalGyrus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset15 <- read.table("Binned_expression_PXD012131_Brain-InferiorParietalLobule.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset16 <- read.table("Binned_expression_PXD012131_Brain-NeoCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset17 <- read.table("Binned_expression_PXD012131_Brain-SuperiorTemporalGyrus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset18 <- read.table("Binned_expression_PXD012131_Brain-Thalamus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset19 <- read.table("Binned_expression_PXD012131_Brain-VisualCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset20 <- read.table("Binned_expression_PXD012755_Brain-CerebellarHemisphericCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset21 <- read.table("Binned_expression_PXD012755_Brain-OccipitalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset22 <- read.table("Binned_expression_PXD015079_Brain-PrefrontalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset23 <- read.table("Binned_expression_syn21443008_Brain_UniPenn-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset24 <- read.table("Binned_expression_syn21444980_Brain_Aging-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset25 <- read.table("Binned_expression_syn3606087_Brain_BLSA-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset26 <- read.table("Binned_expression_syn4624471_Brain_BLSA-Precuneus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset27 <- read.table("Binned_expression_syn6038797_Brain_MountSinai-FrontalPole.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset28 <- read.table("Binned_expression_syn6038852_Brain_ACT-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset29 <- read.table("Binned_expression_syn7204174_Brain_Banner-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset30 <- read.table("Binned_expression_syn7431984_Brain_MayoClinic-TemporalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# Keep the minimum number of genes that are common in all datsets
Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene"), by.y=c("Gene"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(dataset1, dataset2)
merged_data <- Merge_data(merged_data, dataset3)
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
merged_data <- Merge_data(merged_data, dataset26)
merged_data <- Merge_data(merged_data, dataset27)
merged_data <- Merge_data(merged_data, dataset28)
merged_data <- Merge_data(merged_data, dataset29)
merged_data <- Merge_data(merged_data, dataset30)
#colnames(merged_data) <- gsub("AnteriorPitutaryGland", "PitutaryGland", colnames(merged_data), perl=TRUE)
#colnames(merged_data) <- gsub("PitutaryHypophysis", "PitutaryGland", colnames(merged_data), perl=TRUE)

###### HEART

dataset1 <- read.table("Binned_expression_PXD006675_Heart-Aorta.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table("Binned_expression_PXD006675_Heart-AorticValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset3 <- read.table("Binned_expression_PXD006675_Heart-AtrialSeptum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table("Binned_expression_PXD006675_Heart-InferiorVenaCava.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table("Binned_expression_PXD006675_Heart-LeftAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table("Binned_expression_PXD006675_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table("Binned_expression_PXD006675_Heart-MitralValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table("Binned_expression_PXD006675_Heart-PulmonaryArtery.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table("Binned_expression_PXD006675_Heart-PulmonaryValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset10 <- read.table("Binned_expression_PXD006675_Heart-PulmonaryVein.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset11 <- read.table("Binned_expression_PXD006675_Heart-RightAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset12 <- read.table("Binned_expression_PXD006675_Heart-RightVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset13 <- read.table("Binned_expression_PXD006675_Heart-TricuspidValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset14 <- read.table("Binned_expression_PXD006675_Heart-VentricularSeptum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset15 <- read.table("Binned_expression_PXD008722_Heart-LeftAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset16 <- read.table("Binned_expression_PXD008722_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset17 <- read.table("Binned_expression_PXD008722_Heart-RightAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset18 <- read.table("Binned_expression_PXD008934_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset19 <- read.table("Binned_expression_PXD010154_Heart.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# Keep the minimum number of genes that are common in all datsets
Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene"), by.y=c("Gene"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(dataset1, dataset2)
merged_data <- Merge_data(merged_data, dataset3)
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
#merged_data <- Merge_data(merged_data, dataset15)
#merged_data <- Merge_data(merged_data, dataset16)
#merged_data <- Merge_data(merged_data, dataset17)
merged_data <- Merge_data(merged_data, dataset18)
merged_data <- Merge_data(merged_data, dataset19)

filtered_data <- merged_data
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

pca_input <- filtered_data[,-c(1, ncol(filtered_data))]



### Normalisation
norm_input <- filtered_data[,-c(1, ncol(filtered_data))]

sample_names_tissues <- gsub(".*\\.", "", colnames(norm_input), perl=TRUE)
sample_names_tissues <- gsub("breast", "Breast", sample_names_tissues, perl=TRUE)
input_batch_tissues <- data.frame(sample_names_tissues)
input_batch_tissues <- input_batch_tissues %>% group_by(sample_names_tissues) %>% mutate( batchID = cur_group_id() )


sample_names_datasets <- gsub("\\..*", "", colnames(norm_input), perl=TRUE)
sample_names_datasets <- gsub("breast", "Breast", sample_names_datasets, perl=TRUE)
input_batch_datasets <- data.frame(sample_names_datasets)
input_batch_datasets <- input_batch_datasets %>% group_by(sample_names_datasets) %>% mutate( batchID = cur_group_id() )




##### Combat normalisation to remove batch effect
#combat_input <- filtered_data[,-c(1,21,22,33,34,35,40,41,45,46,47, ncol(filtered_data))]
#combat_input <- combat_input[complete.cases(combat_input),]

combat_normalised <- ComBat(norm_input,
                            batch = input_batch_datasets$batchID,
                            mod = NULL,
                            par.prior = TRUE,
                            prior.plots = FALSE,
                            mean.only = TRUE,
                            ref.batch = NULL,
                            BPPARAM = bpparam("SerialParam"))

pca_input <- combat_normalised


####### LIMMA normalisation to remove batch effect
# https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf (page 192)


limma_normalised <-  removeBatchEffect(norm_input, 
                                       batch=input_batch_datasets$batchID, 
                                       batch2=NULL, 
                                       covariates=NULL)



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
  geom_text(aes(label=DatasetID), size=5)+
  labs(x="PC1", y="PC2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #labs(color="Samples") +
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.3,"line"))+
  theme(legend.text=element_text(size=12))+
  guides(col = guide_legend(ncol = 3))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ggtitle("Brain-Samples_bin_values_batch-per-dataset-regions\n[filter: genes detected in at least 50% of samples]\n[number of genes: 3662]")

# IMPORTANT reset key
GeomText$draw_key <- oldK








######## Bin matrix Plot
filtered_data <- filtered_data[1:100, 1:(ncol(filtered_data)-1)]
data_long <- gather(filtered_data, Samples, Bins, colnames(filtered_data)[2]:colnames(filtered_data)[ncol(filtered_data)], factor_key=TRUE)

ggplot(data_long, aes(Samples, Gene)) + 
  geom_tile(aes(fill= Bins), colour="white") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Brain-Samples_binned-by-regions\n[filter: genes detected in at least 50% of all samples]")+
  #scale_fill_viridis(discrete=FALSE) +
  coord_flip()+geom_text(aes(label=Bins), colour="white", size=3)


## Heatmap (code from David)
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=0.4,cex.lab=10)
cordata <- merged_data
cordata[is.na(cordata)] <- 1
my_group1 <- as.numeric(as.factor(gsub("\\..*", "", colnames(merged_data[,-c(1)]), perl=TRUE)))
my_group2 <- as.numeric(as.factor(gsub(".*\\.", "", colnames(merged_data[,-c(1)]), perl=TRUE)))
Datasets <- rainbow(4)[my_group1]
Tissues <- rainbow(15)[my_group2]
colSide <- cbind(Datasets, Tissues)
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(10)
cor_results <- cor(cordata[,-c(1)])
ColSideAnn<-data.frame(Datasets=my_group1,Tissues=my_group2,stringsAsFactors=TRUE)
heatmap3(cor_results, 
         RowSideColors=colSide,
         ColSideColors=colSide,
         col=colMain,
         labRow = FALSE,
         labCol = FALSE,
         ColSideAnn=ColSideAnn,
         ColSideWidth = 10)
legend("right", title="Tissues", legend =levels(as.factor(annot_osw$Disease)),pch=16, pt.cex=1.5, cex=0.75, bty='n',
       col =brewer.pal(9, "Set1")[seq_along(levels(as.factor(annot_osw$Disease)))])
legend("bottomright", title="Type", legend =levels(as.factor(annot_osw$Condition)),pch=16, pt.cex=1.5, cex=0.75, bty='n',
       col =brewer.pal(9, "Set1")[seq_along(levels(as.factor(annot_osw$Condition)))])




#########################
##### Scoring System

# Binning analysis (counting) is done separately for each dataset

# Separate datasets after merging. Separating after merging datasets of same tissue types 
# will provide same number of genes across datasets to compare


#Bins_compare_Median <- function(dataset){
#  median = apply(dataset, 1, median, na.rm = T)
#}
# Getting the median value of all bins for a gene
#dataset1$Median_bins <- Bins_compare_Median(dataset1[,-c(1)])
#dataset2$Median_bins <- Bins_compare_Median(dataset2[,-c(1)])


# Count of Samples per Protein per Bin

# Genes that are not detected in samples after merging 
# (i.e., NA values are also assigned to bin1)
# See above (discussion with Andy Jones)

bin_compare_data <- merged_data
bin_compare_data[is.na(bin_compare_data)] <- 1

Count_samples_per_protein_per_bin <- function(dataset){
  Bin_1 <- apply( dataset, 1, function(x) (length(which(x == 1))/ncol(dataset))*100 )
  Bin_2 <- apply( dataset, 1, function(x) (length(which(x == 2))/ncol(dataset))*100 )
  Bin_3 <- apply( dataset, 1, function(x) (length(which(x == 3))/ncol(dataset))*100 )
  Bin_4 <- apply( dataset, 1, function(x) (length(which(x == 4))/ncol(dataset))*100 )
  Bin_5 <- apply( dataset, 1, function(x) (length(which(x == 5))/ncol(dataset))*100 )
  
  result <- as.data.frame(cbind(Bin_1, Bin_2, Bin_3, Bin_4, Bin_5))
  
  # Getting the bin number which has the highest frequency for a gene
  result <- cbind(result, highestfreqBin = colnames(result)[max.col(result, ties.method="last")])
  
  # Getting the bin number that is equal to or passes the cutoff threshold (50%) for a gene
  result <- cbind(result, cutfoffBin_50percent = colnames(result[,-c(ncol(result))])[apply( result[,-c(ncol(result))], 1, function(x) { ifelse(any(x >= 0.5), which(x >= 0.5), NA ) } )])
}

Bin_counts <- cbind(Gene=bin_compare_data[,c("Gene")], Count_samples_per_protein_per_bin(bin_compare_data[,-c(1)]))
write.table(Bin_counts, "Bin_Density_results_Brain.txt", sep = "\t", row.names = FALSE, quote = FALSE )

Bindata_long <- gather(Bin_counts, Bin, Percent_of_samples, colnames(Bin_counts)[2]:colnames(Bin_counts)[6], factor_key=TRUE)
Bindata_long$Bin <- gsub("Bin_", "", Bindata_long$Bin, perl=TRUE)
Bindata_long$bingroup <- cut(Bindata_long$Percent_of_samples, breaks=c(0,10,20,30,40,50,60,70,80,90,100))

ggplot(Bindata_long[Bindata_long$Percent_of_samples !=0,], aes(x=bingroup, fill=Bin)) + 
  geom_histogram(stat="count")+
  theme_bw()+
  xlab("Percent of samples")+
  ylab("Total counts of occurences of bins for all genes")+
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Sample count per gene per bin - Brain\n(total number of samples = 339)")



#





# Plotting for inspection
plotdata1 <- dataset1[,c("Gene", "Median_bins", "highestfreqBin", "cutfoffBin")]
colnames(plotdata1) <- c("Gene", "PXD006675_Median_bin", "PXD006675_HighestFreqBin", "PXD006675_50percent_CutoffBin")

plotdata2 <- dataset2[,c("Gene", "Median_bins", "highestfreqBin", "cutfoffBin")]
colnames(plotdata2) <- c("Gene", "PXD008722_Median_bin", "PXD008722_HighestFreqBin", "PXD008722_50percent_CutoffBin")

comparion_plotdata <- merge(x=plotdata1, y=plotdata2,
                  by.x=c("Gene"), by.y=("Gene"),
                  all.x=TRUE, all.y=TRUE)

comparion_plotdata <- comparion_plotdata[1:75,]

comparison_plotdata_long <- gather(comparion_plotdata, Samples, Bins, colnames(comparion_plotdata)[2]:colnames(comparion_plotdata)[ncol(comparion_plotdata)], factor_key=TRUE)
comparison_plotdata_long$Bins <- gsub(".*_", "", comparison_plotdata_long$Bins, perl=TRUE)

ggplot(comparison_plotdata_long, aes(factor(Samples, levels = c("Gene",
                                                                "PXD006675_Median_bin",
                                                                "PXD008722_Median_bin",
                                                                "PXD006675_HighestFreqBin",
                                                                "PXD008722_HighestFreqBin",
                                                                "PXD006675_50percent_CutoffBin",
                                                                "PXD008722_50percent_CutoffBin")), Gene)) + 
  geom_tile(aes(fill= Bins), colour="white") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #scale_fill_viridis(discrete=FALSE) +
  xlab("Sample")+
  scale_fill_viridis_d()+
  coord_flip()+geom_text(aes(label=Bins), colour="white", size=3)
