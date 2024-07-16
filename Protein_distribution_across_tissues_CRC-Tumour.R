
######## colorectal cancer Tumour data - Javier ###########

library(ggplot2)
library(dplyr)
library(stats)
library(tidyr)
library(viridis)
library(grid)
library(RColorBrewer)
library(heatmap3)
library(stringr)
library(ggrepel)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
library(sva)

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Javier/CRC_Tumour')


Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene.ID", "Gene.Symbol"), by.y=c("Gene.ID", "Gene.Symbol"), all.x=TRUE, all.y=TRUE)
}

Get_median <- function(dataset){
  dataset$Median = apply(dataset[,-c(1)], 1, median, na.rm = T)
  dataset <- dataset[,c("Gene", "Median"), drop=FALSE]
}


input_data <- read.table(file = "proteinGroups_ppb_final-tissue_names-CRCTumour.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


colnames(input_data) <- gsub("Gene.ID", "GeneID", colnames(input_data))
colnames(input_data) <- gsub("Gene.Symbol", "GeneName", colnames(input_data))
colnames(input_data) <- gsub(".*_", "", colnames(input_data), perl=TRUE)

Gene_info <- input_data[,c("GeneID", "GeneName")]

#Some of the gene entries (ex. IGHA2 has two Ensembl gene ids ENSG00000211890 & ENSG00000276173
#                          ex. IGHV2-70 has two Ensembl gene ids ENSG00000274576 & ENSG00000282453)
# because of this there are duplicate entries of such genes. These are aggregated by taking the median of them

#Merged_input_data <- aggregate(Merged_input_data[,-c(1,2)], list("GeneName"= Merged_input_data$GeneName), median, na.rm =TRUE)


#1. Separate samples by organs and extract genes
Colorectal_adenoma <- input_data[grepl("GeneName|Colorectal.adenoma", colnames(input_data), ignore.case = TRUE)]
Colorectal_tumor <-  input_data[grepl("GeneName|Colorectal.tumor", colnames(input_data), ignore.case = TRUE)]
Colorectal_mucosa <- input_data[grepl("GeneName|Colorectal.mucosa", colnames(input_data), ignore.case = TRUE)]

# Get sample sizes for each organ (i.e., number of MS runs representing each organ)
sample_sizes <- as.data.frame(ncol(Colorectal_adenoma)-1)
sample_sizes <- cbind(sample_sizes, ncol(Colorectal_tumor)-1)
sample_sizes <- cbind(sample_sizes, ncol(Colorectal_mucosa)-1)

colnames(sample_sizes) <- gsub(".*\\(", "", colnames(sample_sizes))
colnames(sample_sizes) <- gsub("\\).*", "", colnames(sample_sizes))
sample_sizes <- as.data.frame(t(sample_sizes))
sample_sizes <- tibble::rownames_to_column(sample_sizes, "Tissues")
colnames(sample_sizes) <- c("Tissues", "Number_of_samples")

# Data to plot distribution of iBAQ values across organs
Colorectal_adenoma_iBAQ_long <- gather(Colorectal_adenoma, Sample, iBAQ, colnames(Colorectal_adenoma)[2]:colnames(Colorectal_adenoma)[ncol(Colorectal_adenoma)], factor_key=TRUE)
Colorectal_tumor_iBAQ_long <- gather(Colorectal_tumor, Sample, iBAQ, colnames(Colorectal_tumor)[2]:colnames(Colorectal_tumor)[ncol(Colorectal_tumor)], factor_key=TRUE)
Colorectal_mucosa_iBAQ_long <- gather(Colorectal_mucosa, Sample, iBAQ, colnames(Colorectal_mucosa)[2]:colnames(Colorectal_mucosa)[ncol(Colorectal_mucosa)], factor_key=TRUE)


All_iBAQ_long <- rbind(Colorectal_adenoma_iBAQ_long,
                       Colorectal_tumor_iBAQ_long,
                       Colorectal_mucosa_iBAQ_long)

All_iBAQ_long <- All_iBAQ_long[complete.cases(All_iBAQ_long),]
All_iBAQ_long$Tissues <- gsub(".*[0-9]\\.", "", All_iBAQ_long$Sample, perl=TRUE)
All_iBAQ_long$Tissues <- gsub("\\.[0-9]", "", All_iBAQ_long$Tissues, perl=TRUE)
All_iBAQ_long$Tissues <- gsub("[0-9]", "", All_iBAQ_long$Tissues, perl=TRUE)
All_iBAQ_long$Tissues <- gsub("\\.", "_", All_iBAQ_long$Tissues, perl=TRUE)
All_iBAQ_long$Tissues <- gsub("CPTAC_", "", All_iBAQ_long$Tissues, perl=TRUE)

All_iBAQ_long <- merge(x=All_iBAQ_long, y=sample_sizes,
                       by.x=c("Tissues"), by.y=c("Tissues"))

All_iBAQ_long$Tissue_samples <- paste(All_iBAQ_long$Tissues, " (", All_iBAQ_long$Number_of_samples, ")", sep="")

All_iBAQ_long$Datasets <- gsub("\\..*", "", All_iBAQ_long$Sample, perl=TRUE)
All_iBAQ_long$Tissue_samples <- gsub("_", " ", All_iBAQ_long$Tissue_samples, perl=TRUE)

ggplot(All_iBAQ_long, aes(x=Tissue_samples, y=iBAQ)) + 
  geom_boxplot() + 
  xlab("Tissues")+
  ylab("protein abundance (ppb)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("iBAQ (CRC Tumor)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  scale_x_discrete(labels = function(x) str_wrap(x, width =15))




# II. Compute the median from all samples of each Tissue separately
Tissue_median <- as.data.frame(cbind(GeneName = Colorectal_adenoma$GeneName, Colorectal_adenoma = apply(as.data.frame(Colorectal_adenoma[,-c(1)]), 1, median, na.rm = T)))
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Colorectal_tumor$GeneName, Colorectal_tumor = apply(Colorectal_tumor[,-c(1)], 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Tissue_median <- merge(x=Tissue_median, y=cbind(GeneName = Colorectal_mucosa$GeneName, Colorectal_mucosa = apply(as.data.frame(Colorectal_mucosa[,-c(1)]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)

write.table(Tissue_median, file = paste("Gene_distribution_in_Tissues-GeneNames-Median_intensities-CRC_Tumor.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

Count <- data.frame(Present_in_number_of_samples = apply(Tissue_median[2:ncol(Tissue_median)], 1, function(x) length(which(x != "NaN")) ))
#Count <- Count[Count$Present_in_number_of_samples !=0,,drop=FALSE]
Count$Sample_percentage <- ((Count$Present_in_number_of_samples)/(ncol(Tissue_median)-1)) *100
Count <- group_by(Count, Present_in_number_of_samples) %>% mutate(Total_number_of_sample_occurences =n()) %>% mutate(Gene_percent = (Total_number_of_sample_occurences / nrow(Count))*100)

Tissue_median <- cbind(Tissue_median, Count)

Tissue_median <- merge(x=Gene_info, y=Tissue_median,
                       by.x=c("GeneName"), by.y=c("GeneName"),
                       all.x=FALSE, all.y=FALSE)

plotdata <- unique(Count)


ggplot(plotdata[plotdata$Present_in_number_of_samples !=0,], aes(x=Present_in_number_of_samples, y=Gene_percent)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(limits=c("1","2","3"))+
  xlab("Present in number of tissues")+
  ylab("% of identified 'canonical proteins'")+
  theme_bw()+
  ggtitle("Distribution of 'canonical proteins' across tissues\n - (CRC Tumor)")

ggplot(plotdata[plotdata$Present_in_number_of_samples !=0,], aes(x=Present_in_number_of_samples, y=Total_number_of_sample_occurences)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(limits=c("1","2","3"))+
  xlab("Present in number of tissues")+
  ylab("Number of canonical proteins")+
  theme_bw()+
  ggtitle("Distribution of 'canonical proteins' across tissues\n (CRC Tumor)")


###### Gene distribution per organ
gene_counts_per_tissue <- as.data.frame(sapply(Tissue_median[,-c(1:2)], function(x) sum(!is.na(x))))
colnames(gene_counts_per_tissue) <- "Number_of_identified_genes"
gene_counts_per_tissue <- tibble::rownames_to_column(gene_counts_per_tissue, "Tissues")
gene_counts_per_tissue$Tissues <- gsub("_median", "", gene_counts_per_tissue$Tissues, perl=TRUE)
gene_counts_per_tissue$Percentage <- (gene_counts_per_tissue$Number_of_identified_genes/nrow(Tissue_median))*100

gene_counts_per_tissue <- merge(x=gene_counts_per_tissue, y=sample_sizes,
                                by.x=c("Tissues"), by.y="Tissues")


gene_counts_per_tissue$Tissue_samples <- paste(gene_counts_per_tissue$Tissue, " (", gene_counts_per_tissue$Number_of_samples, ")", sep="")
gene_counts_per_tissue$Tissue_samples <- gsub("_", " ", gene_counts_per_tissue$Tissue_samples, perl=TRUE)

write.table(gene_counts_per_tissue, file = paste("Gene_distribution_in_tissues_plot-CRC_Tumor.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

ggplot(gene_counts_per_tissue, aes(x=Tissue_samples, y=Percentage)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("% identified canonical proteins")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across tissues\n (CRC-Tumor)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  scale_x_discrete(labels = function(x) str_wrap(x, width =15))

ggplot(gene_counts_per_tissue, aes(x=Tissue_samples, y=Number_of_identified_genes)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across tissues\n (CRC-Tumor)")+
  theme(axis.text=element_text(size=12,),
        axis.title=element_text(size=12))+
  scale_x_discrete(labels = function(x) str_wrap(x, width =15))

##### Gene distribution per dataset
All_datasets <- input_data

dataset_sample_names <- colnames(All_datasets[-c(1,2)])
colnames(All_datasets) <- gsub("\\..*", "", colnames(All_datasets), perl=TRUE)


Datasets_median <- as.data.frame(cbind(GeneName = All_datasets$GeneName, CPTAC_median = apply(as.data.frame(All_datasets[,grep("CPTAC", colnames(All_datasets))]), 1, median, na.rm = T)))
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD001676_median = apply(as.data.frame(All_datasets[,grep("PXD001676", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD002137_median = apply(as.data.frame(All_datasets[,grep("PXD002137", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD014511_median = apply(as.data.frame(All_datasets[,grep("PXD014511", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)
Datasets_median <- merge(x=Datasets_median, y=cbind(GeneName = All_datasets$GeneName, PXD019504_median = apply(as.data.frame(All_datasets[,grep("PXD019504", colnames(All_datasets))]), 1, median, na.rm = T)), by.x="GeneName", by.y="GeneName", all.x=TRUE, all.y=TRUE)


# Gene distribution per dataset
gene_counts_per_dataset <- as.data.frame(sapply(Datasets_median[,-c(1)], function(x) sum(!is.na(x))))
colnames(gene_counts_per_dataset) <- "Number_of_identified_genes"
gene_counts_per_dataset <- tibble::rownames_to_column(gene_counts_per_dataset, "Datasets")
gene_counts_per_dataset$Datasets <- gsub("_median", "", gene_counts_per_dataset$Datasets, perl=TRUE)
gene_counts_per_dataset$Percentage <- (gene_counts_per_dataset$Number_of_identified_genes/nrow(Datasets_median))*100

# Number of tissues per dataset
dataset_tissues <- data.frame(Datasets= "CPTAC", Tissue_count = length(unique(gsub("CPTAC.", "", dataset_sample_names[grep("CPTAC", dataset_sample_names)], perl=TRUE))) )
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD001676", Tissue_count = length(unique(gsub("PXD[0-9].", "", dataset_sample_names[grep("PXD001676", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD002137", Tissue_count = length(unique(gsub("PXD[0-9].", "", dataset_sample_names[grep("PXD002137", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD014511", Tissue_count = length(unique(gsub("PXD[0-9].", "", dataset_sample_names[grep("PXD014511", dataset_sample_names)], perl=TRUE))) ))
dataset_tissues <- rbind(dataset_tissues, data.frame(Datasets= "PXD019504", Tissue_count = length(unique(gsub("PXD[0-9].", "", dataset_sample_names[grep("PXD019504", dataset_sample_names)], perl=TRUE))) ))



gene_counts_per_dataset <- merge(x=gene_counts_per_dataset, y=dataset_tissues,
                                 by.x=c("Datasets"), by.y=c("Datasets"), all.x=FALSE, all.y=FALSE)

gene_counts_per_dataset$Datasets_tissues <- paste(gene_counts_per_dataset$Datasets, " (", gene_counts_per_dataset$Tissue_count, ")", sep="")

write.table(gene_counts_per_dataset, file = paste("Gene_distribution_in_datasets_plot-CRC_Tumor.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )


ggplot(gene_counts_per_dataset, aes(x=Datasets_tissues, y=Number_of_identified_genes)) + 
  geom_bar(stat="identity") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across datasets\n (CRC-Tumor)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

# Plot distribution of iBAQ values across datasets
All_iBAQ_long <- merge(x=All_iBAQ_long, y=dataset_tissues,
                       by.x=c("Datasets"), by.y=c("Datasets"))

All_iBAQ_long$Dataset_tissues <- paste(All_iBAQ_long$Datasets, " (", All_iBAQ_long$Tissue_count, ")", sep="")


ggplot(All_iBAQ_long, aes(x=Dataset_tissues, y=iBAQ)) + 
  geom_boxplot() + 
  xlab("Datasets")+
  ylab("protein abundance (ppb)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("iBAQ (CRC-Tumor)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))


# Relation between sample size and gene coverage
scatterplot <- gene_counts_per_tissue[,c("Percentage", "Tissues")]
scatterplot$Name <- gsub(" .*", "", scatterplot$Tissues, perl=TRUE)
scatterplot$Tissues <- gsub(".*\\(|\\)", "", scatterplot$Tissues, perl=TRUE)


ggplot(scatterplot, aes(x=as.numeric(Tissues), y=Percentage)) + 
  geom_point() + 
  geom_smooth(method='lm')+
  xlab("sample size")+
  ylab("% canonical protein coverage")+
  theme_bw()+
  scale_x_log10()+
  geom_label_repel(aes(label = Name), size = 3, max.overlaps = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Relation between sample size and 'canonical protein' coverage")+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))


###################################################################
############################## E N D ##############################
###################################################################
##### don't use the PCA here, use the method in rank-binning script  



Brain_long <- gather(Brain, Samples, Intensities, colnames(Brain[2]):colnames(Brain[ncol(Brain)]), factor_key=TRUE)
Brain_long$Tissues <- gsub(".*\\.", "", Brain_long$Samples, perl=TRUE)

# If gene is identified in at least one tissue sample then that gene is counted as identified
Brain_long <- Brain_long %>% group_by(GeneName, Tissues) %>% mutate(is_gene_NA_in_all_tissues = all(is.na(Intensities)))

filtered_data_long <- Brain_long[Brain_long$is_gene_NA_in_all_tissues == "FALSE",]



#Heart_long <- gather(Heart, Samples, Intensities, colnames(Heart[2]):colnames(Heart[ncol(Heart)]), factor_key=TRUE)
#Heart_long$Tissues <- gsub(".*\\.", "", Heart_long$Samples, perl=TRUE)

# If gene is identified in at least one tissue sample then that gene is counted as identified
#Heart_long <- Heart_long %>% group_by(GeneName, Tissues) %>% mutate(is_gene_NA_in_all_tissues = all(is.na(Intensities)))

#filtered_data_long <- Heart_long[Heart_long$is_gene_NA_in_all_tissues == "FALSE",]




### 
# The undetected genes are also put into bin1 
# (After discussing with Andy Jones [20/01/2021]: since undetected genes could be below detection threshold due to their low expression or abundance)
#replace missing values (NA) of those genes that were detected in at least one tissue sample
#filtered_data_long[is.na(filtered_data_long)] <- 0

filtered_data <- spread(filtered_data_long[,-c(4,5)], Samples, Intensities)

# count number of samples for each gene that have non NA values (i.e., detected)
non_missing_sample_count_percentage <- data.frame(non_missing_sample_count_percentage= apply(filtered_data[2:ncol(filtered_data)], 1, function(x) (length(which(!is.na(x))))/(ncol(filtered_data)-1)*100 ))

filtered_data <- cbind(filtered_data, non_missing_sample_count_percentage)

# Filter genes that are present in at least 33%, 50% or 75% of samples
filtered_data <- filtered_data[filtered_data$non_missing_sample_count_percentage >= 50,]


# PCA of only intensity values without batch correction
pca_input <- filtered_data[,-c(1,ncol(filtered_data))]
pca_input[is.na(pca_input)] <- 0



#### Normalisation 
norm_input <- filtered_data[,-c(1,ncol(filtered_data))]
norm_input[is.na(norm_input)] <- 0

# Remove tissue samples with only 1 replicate
# norm_input <- filtered_data[,-c(1,21,22,33,34,35,40,41,45,46,47, ncol(filtered_data))]

sample_names_tissues <- gsub(".*\\.", "", colnames(norm_input), perl=TRUE)
sample_names_tissues <- gsub("breast", "Breast", sample_names_tissues, perl=TRUE)
input_batch_tissues <- data.frame(sample_names_tissues)
input_batch_tissues <- input_batch_tissues %>% group_by(sample_names_tissues) %>% mutate( batchID = cur_group_id() )


sample_names_datasets <- gsub("\\..*", "", colnames(norm_input), perl=TRUE)
sample_names_datasets <- gsub("breast", "Breast", sample_names_datasets, perl=TRUE)
input_batch_datasets <- data.frame(sample_names_datasets)
input_batch_datasets <- input_batch_datasets %>% group_by(sample_names_datasets) %>% mutate( batchID = cur_group_id() )


#### Combat normalisation to remove batch effects
combat_normalised <- ComBat(as.matrix(norm_input),
                            batch = input_batch_datasets$batchID,
                            mod = NULL,
                            par.prior = TRUE,
                            prior.plots = FALSE,
                            mean.only = TRUE,
                            ref.batch = NULL,
                            BPPARAM = bpparam("SerialParam"))

# PCA of combat normalised intensity values
pca_input <- combat_normalised


####### LIMMA normalisation to remove batch effect
# https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf (page 192)


limma_normalised <-  removeBatchEffect(norm_input, 
                                       batch=input_batch_datasets$batchID, 
                                       batch2=NULL, 
                                       covariates=NULL)

# PCA of limma normalised intensity values
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
  geom_text(aes(label=DatasetID))+
  labs(x="PC1", y="PC2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #labs(color="Samples") +
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.3,"line"))+
  guides(col = guide_legend(ncol = 3))+
  ggtitle("Heart-Samples_binned-by-regions batch-per-datasets\n[filter: genes detected in at least 50% of samples]")

# IMPORTANT reset key
GeomText$draw_key <- oldK


#foo <- Merged_input_data[,c(1:6)] %>%
#       group_by(GeneID, GeneName) %>% 
#       mutate(PXD000547 = all(is.na(c(PXD000547.Sample1.CorpusCallosum, PXD000547.Sample2.CorpusCallosum))), 
#              PXD000548 = all(is.na(c(PXD000548.Sample1.AnteriorTemporalLobe, PXD000548.Sample2.AnteriorTemporalLobe))))
