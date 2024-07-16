
combine_datsets <- function(dataset1, dataset2){
  combined <- merge(x=dataset1, y=dataset2, 
                    by.x=c("Gene.ID", "Gene.Name"), by.y=c("Gene.ID", "Gene.Name"), 
                    all.x=TRUE, all.y=TRUE)
}

#Already ppb (FOT) normalised intensity values

tmp  <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD010271/proteinGroups_ppb_final-tissue_names.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
tmp1 <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD008934/proteinGroups_ppb_final-tissue_names.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
tmp2 <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD001325/proteinGroups_ppb_final-tissue_names.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
tmp3 <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD012431/proteinGroups_ppb_final-tissue_names.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

joined_dataset <- combine_datsets(tmp, tmp1)
joined_dataset <- combine_datsets(joined_dataset, tmp2)
joined_dataset <- combine_datsets(joined_dataset, tmp3)


joined_dataset[is.na(joined_dataset)] <- 0
colnames(joined_dataset) <- gsub(".*_", "", colnames(joined_dataset), perl=TRUE)



all_data_Pancreas <- joined_dataset[ ,c(1:2, grep("Pancreas", colnames(joined_dataset),ignore.case=TRUE)) ]
all_data_Liver <- joined_dataset[ ,c(1:2, grep("Liver", colnames(joined_dataset), ignore.case=TRUE)) ]
all_data_Ovary <- joined_dataset[ ,c(1:2, grep("Ovary", colnames(joined_dataset), ignore.case=TRUE)) ]
all_data_Substantia_nigra <- joined_dataset[ ,c(1:2, grep("Substantia.Nigra", colnames(joined_dataset), ignore.case=TRUE)) ]
all_data_Heart <- joined_dataset[ ,c(1:2, grep("Heart", colnames(joined_dataset), ignore.case=TRUE)) ]
all_data_Breast <- joined_dataset[ ,c(1:2, grep("Breast", colnames(joined_dataset), ignore.case=TRUE)) ]

gene_names <- joined_dataset[,c(2)] 



convert_to_scaled_data <- function(data, bins) {
  apply(data, 1, function(sample) {
    ceiling((sample / max(sample, na.rm = T)) * bins)
  })
}

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
    
    # Do some stuff.
  })
}

Binned_data_Pancreas <- Binning(all_data_Pancreas[-c(1:2)], gene_names)

Bins_5_Pancreas <- Binned_data_Pancreas$bins_5
Bins_10_Pancreas<- Binned_data_Pancreas$bins_10

Binned_data_Breast <- Binning(all_data_Breast[-c(1:2)], gene_names)

Bins_5_Breast <- Binned_data_Breast$bins_5
Bins_10_Breast<- Binned_data_Breast$bins_10

Median_Bins_5_Breast <- apply(Bins_5_Breast,2,median)
Median_Bins_5_Breast <- data.frame(Gene_ID=names(Median_Bins_5_Breast), Median_bins=unlist(Median_Bins_5_Breast))
Median_Bins_5_Breast$Tissue <- rep("Breast",nrow(Median_Bins_5_Breast))

