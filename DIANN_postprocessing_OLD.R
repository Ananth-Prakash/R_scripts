# DIANN postprocessing
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggVennDiagram)
library(viridis)


setwd('/Users/ananth/Documents/DIANN/')

#Read input file
#"PXD012254", #control samples only
#"PXD001506", #control samples only
#"PXD001764", #control samples only
#"PXD025705", #control samples only
#"PXD004684", #control samples only
#"PXD032076", #control samples only
#"PXD017051", #cases only
#"PXD025431", #control samples only
#"PXD002732", #control samples only
#"PXD000672", #control samples only
#"PXD004873", #control samples only
#"PXD019594", #control samples only
#"PXD022872", #control samples only
#"PXD031419", #control samples only
#"PXD034908", #DIA samples only
#"PXD018430", #all samples, authors contacted, dataset dropped
#"PXD027906", #all samples, authors contacted, dataset dropped
#"PXD033060", #control samples only
#"PXD018830", #control samples only
#"PXD039665", #control samples only
#"PXD022872", #control samples only
#"PXD022952", #authors contacted
#"PXD018678", #authors contacted
#"PXD028854", #authors contacted
#"PXD014943", #authors contacted
#"PXD029359", #authors contacted

DIA_datasets <- data.frame(Dataset = c(
 "PXD012254",
 "PXD001506",
 "PXD001764",
 "PXD025705",
 "PXD004684",
 "PXD032076",
 "PXD025431",
 "PXD002732",
 "PXD000672",
 "PXD004873",
 "PXD019594",
 "PXD022872",
 "PXD034908",
 "PXD018830",
 "PXD039665",
 "PXD033060",
 "PXD031419"
))

#Rootdata <- "SwissProt"
Rootdata <- "RefProt"

combined_reports <- list()


for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  print(datasetID)
  tmp  <- read.table(file=paste(datasetID, Rootdata, "report_controlsamples.tsv", sep="/"), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  tmp$File.Name <- gsub(".*PXD","PXD",tmp$File.Name,perl=TRUE)
  tmp$File.Name <- gsub("/submitted.*","",tmp$File.Name,perl=TRUE)
  
  # don't consider peptides which map to more than one gene (gene groups)
  # peptide mappings should be unique to a gene
  tmp <- tmp[!grepl(";",tmp$Genes),]
  tmp <- tmp[tmp$Genes!= "",]
  #Remove contaminants
  tmp <- tmp[!grepl("SWISS-PROT",tmp$Genes),]
  
  tmp <- tmp[,c("File.Name","Run","Protein.Group","Protein.Names","Genes","Genes.MaxLFQ","Genes.Normalised","Stripped.Sequence")]
  tmp$Genes.Normalised <- as.numeric(tmp$Genes.Normalised)
  tmp$Genes.MaxLFQ <- as.numeric(tmp$Genes.MaxLFQ)

  combined_reports[[i]] <- tmp
}

reports_submatrix <- bind_rows(combined_reports)
reports_submatrix <- reports_submatrix[reports_submatrix$Genes!="" & reports_submatrix$Genes!="(Bos" & reports_submatrix$Genes!="(S.avidinii)", ]

peptide_counts_per_gene <- reports_submatrix  %>%
  group_by(File.Name,Run,Protein.Group,Protein.Names,Genes) %>%  #Protein.Group can be removed to give accurate count
  summarise(Total_peptides = n(),
            Unique_peptides = n_distinct(Stripped.Sequence))

peptide_counts_per_gene <- peptide_counts_per_gene[peptide_counts_per_gene$File.Name != "",]

#ggplot(peptide_counts_per_gene[peptide_counts_per_gene$File.Name=="PXD000672",], aes(x=Total_peptides, colour=Run, after_stat(count)))+
ggplot(peptide_counts_per_gene, aes(x=Total_peptides, colour=Run, after_stat(count)))+
  geom_density()+
  scale_x_log10(breaks = c(1,2,5,10,25,100,1000))+
  xlab("Total peptidoforms per gene")+
  theme_bw()+
  ggtitle("File: report.tsv, RefProt(100K set)\nNumber of peptidoforms (all) per gene")+
  theme(legend.position="none")+
  facet_wrap(~File.Name, scales="free_y")
  
ggplot(peptide_counts_per_gene, aes(x=Unique_peptides, colour=Run, after_stat(count)))+
  geom_density()+
  scale_x_log10(breaks = c(1,2,5,10,25,100,1000))+
  xlab("Unique peptidoforms per gene")+
  theme_bw()+
  ggtitle("File: report.tsv, RefProt(100K set)\nNumber of peptidoforms (unique) per gene")+
  theme(legend.position="none")+
  facet_wrap(~File.Name, scales="free_y")



# Check if to filter by "Total peptides" or "Unique peptides"?
# Total_peptides  = Total number of peptidoforms mapped to a gene (can have duplicates)
# Unique_peptides = Total number of unique peptidoforms only that are mapped to a gene (does not have duplicates) 

More_than_two_peptides <-  peptide_counts_per_gene[peptide_counts_per_gene$Total_peptides > 1, ]
#More_than_two_peptides <- peptide_counts_per_gene[peptide_counts_per_gene$Unique_peptides > 1, ]

rm(combined_reports)

all_genes <- data.frame(More_than_two_peptides$Genes)
genes_2ormore_pept <- data.frame(Genes=unique(all_genes$More_than_two_peptides.Genes))


all_proteins <- More_than_two_peptides[,c("Genes","Protein.Group","Protein.Names")]
proteins_2ormore_pept <- unique(all_proteins)

DIA_tissues_dataset_annot <- read.table(file="DIA_tissues_datasets.tsv", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)


#################################

#1 read unique_genes_matrix files
combined_GeneLFQs_2ormore_pept <- genes_2ormore_pept 
for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  print(datasetID)
  tmp  <- read.table(file=paste(datasetID, Rootdata, "report.unique_genes_matrix_controlsamples.tsv", sep="/"), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", check.names=FALSE)
  colnames(tmp) <- gsub(".*\\/","",colnames(tmp),perl=TRUE)
  colnames(tmp) <- paste(datasetID,colnames(tmp), sep="~")
  colnames(tmp)[1] <- c("Genes")
  # Remove any mappings to more than one gene
  tmp <- tmp[!grepl(";",tmp$Genes),]
  tmp <- tmp[tmp$Genes!= "",]
  #Remove contaminants
  tmp <- tmp[!grepl("SWISS-PROT",tmp$Genes),]
  
  ### merge tables of all datasets and then create a long format outside
  combined_GeneLFQs_2ormore_pept <- merge(x=combined_GeneLFQs_2ormore_pept, y=tmp,
                                by.x=c("Genes"), by.y=c("Genes"),
                                all.x=TRUE, all.y=FALSE)
}

combined_GeneLFQs_2ormore_pept_Long <- gather(combined_GeneLFQs_2ormore_pept, Sample, LFQ, colnames(combined_GeneLFQs_2ormore_pept)[2]:colnames(combined_GeneLFQs_2ormore_pept)[ncol(combined_GeneLFQs_2ormore_pept)], factor_key=TRUE)
combined_GeneLFQs_2ormore_pept_Long$Dataset <- gsub("~.*","",combined_GeneLFQs_2ormore_pept_Long$Sample, perl=TRUE)

ggplot(combined_GeneLFQs_2ormore_pept_Long, aes(x=LFQ, color=Sample))+
  scale_x_log10()+
  geom_density()+
  theme_bw()+
  ggtitle("File: report.unique_genes_matrix.tsv, RefProt(100K set)\nFilter: Total number of peptides per gene > 1")+
  theme(legend.position="none")+
  facet_wrap(~Dataset, scales = "free")

# QC plot. PXD017051 has outlier samples
# PXD017051~16-T(N)105.mzML & 
# PXD017051~17-T(N)072.mzML

ggplot(combined_GeneLFQs_2ormore_pept_Long, aes(x=Sample, y= LFQ, color=Sample))+
#ggplot(combined_GeneLFQs_2ormore_pept_Long[combined_GeneLFQs_2ormore_pept_Long$Dataset == "PXD017051",], aes(x=Sample, y= LFQ, color=Sample))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ggtitle("File: report.unique_genes_matrix.tsv, RefProt(100K set)\nFilter: Total number of peptides per gene > 1")+
  theme(legend.position="none")+
  facet_wrap(~Dataset, scales = "free")

# Remove outlier samples
combined_GeneLFQs_2ormore_pept_Long <- combined_GeneLFQs_2ormore_pept_Long[combined_GeneLFQs_2ormore_pept_Long$Sample != "PXD017051~16-T(N)105.mzML" & combined_GeneLFQs_2ormore_pept_Long$Sample != "PXD017051~17-T(N)072.mzML",]

combined_GeneLFQs_2ormore_pept_Long <- merge(x=DIA_tissues_dataset_annot, y=combined_GeneLFQs_2ormore_pept_Long,
                                         by.x=c("Dataset"), by.y=c("Dataset"),
                                         all.x=TRUE, all.y=FALSE)

#################################
#2 read protein groups matrix files
combined_ProteinLFQs_2ormore_pept <- proteins_2ormore_pept

for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  print(datasetID)
  tmp <- read.table(file=paste(datasetID, Rootdata, "report.pg_matrix_controlsamples.tsv", sep="/"), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", check.names=FALSE)
  # Remove protein group column
  tmp <- subset(tmp, select = -c(Protein.Group, First.Protein.Description,Protein.Ids,Protein.Names))

  # Remove any mappings to more than one gene ID
  tmp <- tmp[!grepl(";",tmp$Genes),]
  tmp <- tmp[tmp$Genes!= "",]
  #Remove contaminants
  tmp <- tmp[!grepl("SWISS-PROT",tmp$Genes),]
  
  # Some genes ids have duplicated protein groups, therefore removing duplicate entries by averaging (median) LFQs on gene ids
  tmp <- aggregate(tmp[ , 2:ncol(tmp)], list("Genes" = tmp$Genes), median, na.rm =TRUE)
  
  colnames(tmp) <- gsub(".*\\/","",colnames(tmp),perl=TRUE)
  colnames(tmp) <- paste(datasetID,colnames(tmp), sep="~")
  colnames(tmp) <- gsub(".*~Genes","Genes",colnames(tmp), perl=TRUE)
  colnames(tmp) <- gsub(".*~Protein.Ids","Protein.Ids",colnames(tmp), perl=TRUE)
  colnames(tmp) <- gsub(".*~Protein.Names","Protein.Names",colnames(tmp), perl=TRUE)
  
  print(colnames(tmp))
  
  ### merge tables of all datasets and then create a long format outside
  combined_ProteinLFQs_2ormore_pept <- merge(x=combined_ProteinLFQs_2ormore_pept, y=tmp,
                                          by.x=c("Genes"), by.y=c("Genes"),
                                          all.x=TRUE, all.y=FALSE)
}
combined_ProteinLFQs_2ormore_pept_Long <- gather(combined_ProteinLFQs_2ormore_pept, Sample, LFQ, colnames(combined_ProteinLFQs_2ormore_pept)[4]:colnames(combined_ProteinLFQs_2ormore_pept)[ncol(combined_ProteinLFQs_2ormore_pept)], factor_key=TRUE)
combined_ProteinLFQs_2ormore_pept_Long$Dataset <- gsub("~.*","",combined_ProteinLFQs_2ormore_pept_Long$Sample, perl=TRUE)
#Remove outlier samples
combined_ProteinLFQs_2ormore_pept_Long <- combined_ProteinLFQs_2ormore_pept_Long[combined_ProteinLFQs_2ormore_pept_Long$Sample != "PXD017051~16-T(N)105.mzML" & combined_ProteinLFQs_2ormore_pept_Long$Sample != "PXD017051~17-T(N)072.mzML",]

combined_ProteinLFQs_2ormore_pept_Long <- combined_ProteinLFQs_2ormore_pept_Long[ grep(";", combined_ProteinLFQs_2ormore_pept_Long$Protein.Names, invert = TRUE) , ]

combined_ProteinLFQs_2ormore_pept_Long <- unique(combined_ProteinLFQs_2ormore_pept_Long %>%
                                                   separate_longer_delim(c(Genes, Protein.Group, Protein.Names), delim = ";"))

combined_ProteinLFQs_2ormore_pept_Long <- merge(x=DIA_tissues_dataset_annot, y=combined_ProteinLFQs_2ormore_pept_Long,
                                             by.x=c("Dataset"), by.y=c("Dataset"),
                                             all.x=FALSE, all.y=FALSE)


ggplot(combined_ProteinLFQs_2ormore_pept_Long, aes(x=LFQ, color=Sample))+
  scale_x_log10()+
  geom_density()+
  theme_bw()+
  ggtitle("File: report.pg_matrix_controlsamples.tsv, RefProt(100K set)\nFilter: Total number of peptides per gene > 1")+
  theme(legend.position="none")+
  facet_wrap(~Dataset, scales = "free")

#####################################################################
# compare results with DDA for specific tissues/organs
# Tissues       DIA                   DDA
# Colon       PXD012254               PXD010154, PXD001608, PXD002029
# Duodenum    PXD001506               PXD010154
# Heart       PXD019594               PXD010154, PXD006675, PXD008934
# Liver       PXD025705, PXD004873    PXD010154, PXD010271
# Lung        PXD004684               PXD010154
# Pancreas    PXD032076, PXD017051    PXD010154, PXD010271
# Thyroid     PXD002732               PXD010154
# Uterus      PXD027906               PXD010154
# Kidney      PXD000672               PXD010154
# Esophagus   PXD001764               PXD010154
# Brain       PXD033060, PXD022872    PXD010154, PXD005819, PXD004143, PXD006233, .....
#####################################################################


### Read DDA iBAQ values
DDA_all <- read.table(file="/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/JPR/SupplementaryFiles/SupplementaryTable_2.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
Protein_geneid_map <- DDA_all[,c("Majority.protein.IDs","Gene.Symbol","EnsemblID")]

DDA_all <- DDA_all[,c("Gene.Symbol","Majority.protein.IDs", "EnsemblID","Colon","Duodenum","Esophagus","Heart","Kidney","Liver","Lung","Pancreas","Thyroid","UterineEndometrium","Brain")]

DDA_all <- DDA_all[!grepl(";", DDA_all$EnsemblID),]
DDA_all <- DDA_all[!grepl(";", DDA_all$Gene.Symbol),]

DDA_all_separate <- unique(DDA_all %>%
                    separate_longer_delim(c(Majority.protein.IDs, Gene.Symbol, EnsemblID), delim = ";"))

DDA_all_aggregate <- aggregate(DDA_all_separate[ ,4:ncol(DDA_all_separate)], list("Gene.ID" = DDA_all_separate$EnsemblID, "Gene.Symbol" = DDA_all_separate$Gene.Symbol, "Protein.Ids" = DDA_all_separate$Majority.protein.IDs), median, na.rm =TRUE)

DDA_all_long <- gather(DDA_all_aggregate, Samples, ppb.iBAQ, colnames(DDA_all_aggregate)[4]:colnames(DDA_all_aggregate)[ncol(DDA_all_aggregate)])

DDA_all_long$Tissues <- DDA_all_long$Samples
DDA_all_long$Tissues <- gsub("UterineEndometrium","Uterus",DDA_all_long$Tissues)

colnames(DDA_all_long)[5] <- "Median_ppb.iBAQ"

#Median LFQs and iBAQs (over samples) across tissues and then compare between DDA and DIA
#Should not compare sample vs sample

# Already carried out in earlier steps
# DDA_all_median_tissues <- DDA_all_long  %>%
#  group_by(Gene.Symbol,Gene.ID,Tissues,Protein.Ids) %>%
#  summarise(Median_ppb.iBAQ = median(Median_ppb.iBAQ, na.rm=TRUE))


combined_ProteinLFQs_2ormore_pept_median_tissues <- combined_ProteinLFQs_2ormore_pept_Long %>%
  group_by(Tissues, Genes, Protein.Group) %>%
  summarise(Median_LFQ = median(LFQ, na.rm=TRUE))

combined_ProteinLFQs_2ormore_pept_median_tissues$Tissues <- gsub("Esophageal epithelium","Esophagus", combined_ProteinLFQs_2ormore_pept_median_tissues$Tissues)

merged_DDA_DIA <- merge(x=combined_ProteinLFQs_2ormore_pept_median_tissues, y=DDA_all_long,
                        by.x=c("Genes","Tissues","Protein.Group"), by.y=c("Gene.Symbol","Tissues","Protein.Ids"),
                        all.x=FALSE, all.y=FALSE)


ggplot(merged_DDA_DIA, aes(x=Median_LFQ, y=Median_ppb.iBAQ)) + 
  geom_point(size=0.5, alpha=0.2) + 
  scale_x_log10()+
  scale_y_log10()+
  xlab("Median LFQ (DIA)")+
  ylab("Median ppb.iBAQ (DDA)")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x, size=0.4) +
  ggtitle("Protein abundance comparison\nppb.iBAQ(DDA-RefProt) v/s LFQ(DIA-RefProt 100KSet)")+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=12))+
  facet_wrap(~Tissues, scales = "free")

#Compare genes identified in DDA and DIA
DDA_genes <- data.frame(Gene.Symbol=unique(DDA_all_long[complete.cases(DDA_all_long),c("Gene.Symbol"), drop=FALSE]))
DIA_genes <- combined_ProteinLFQs_2ormore_pept[,c("Genes"), drop=FALSE]

VennPlot <- list(DDA_genes=DDA_genes$Gene.Symbol, DIA_genes=DIA_genes$Genes)
ggVennDiagram(VennPlot)+
  theme(legend.position = "none")+
  ggtitle("Genes identified in DDA-RefProt v/s DIA-RefProt 100K set\nin common organs\n")



##############################################################################
## LFQ values from DIA-NN converted to iBAQ ppb using iBAQpy
#  Correlation with iBAQpy
##############################################################################
# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")
#BiocManager::install("GenomicFeatures")

library(mygene)

ibaqpy_ppb <- list()

for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  print(datasetID)
  foo  <- read.table(file=paste(datasetID, "/RefProt/", datasetID,"-ibaq.tsv", sep=""), quote = "\"", header = TRUE, sep = ",", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  foo$File.Name <- rep(datasetID,nrow(foo))

  
  foo$Ibaq <- as.numeric(foo$Ibaq)
  foo$IbaqNorm <- as.numeric(foo$IbaqNorm)
  foo$IbaqLog <- as.numeric(foo$IbaqLog)
  foo$IbaqPpb <- as.numeric(foo$IbaqPpb)
  
  ibaqpy_ppb[[i]] <- foo
}

ibaqpy_ppb <- bind_rows(ibaqpy_ppb)

Protein_geneid_map_one_to_one <- Protein_geneid_map[ grep(";", Protein_geneid_map$EnsemblID, invert = TRUE) , ]

Protein_geneid_map_one_to_one <- unique(Protein_geneid_map_one_to_one %>%
  separate_longer_delim(c(Majority.protein.IDs, Gene.Symbol, EnsemblID), delim = ";"))

ibaqpy_ppb <- merge(x=ibaqpy_ppb, y=Protein_geneid_map_one_to_one,
                             by.x=c("ProteinName"), by.y=c("Majority.protein.IDs"),
                             all.x=FALSE, all.y=FALSE)

# Average intensities over geneid and dataset. 
# Each DIA dataset here represents only one tissue but multiple samples with normal(healthy) condition
# therefore averaging over all samples and condition
ibaqpy_ppb <- subset(ibaqpy_ppb, select = -c(ProteinName,SampleID,Condition) )
ibaqpy_ppb_aggregate <- ibaqpy_ppb %>%
                                  group_by(File.Name, Gene.Symbol, EnsemblID) %>%
                                  summarise(median_Ibaq = median(Ibaq), median_IbaqNorm = median(IbaqNorm),
                                            median_IbaqLog = median(IbaqLog), median_IbaqPpb = median(IbaqPpb))
  
ibaqpy_ppb_aggregate <- merge(x=DIA_tissues_dataset_annot, y=ibaqpy_ppb_aggregate,
                                       by.x=c("Dataset"), by.y=c("File.Name"),
                                       all.x=FALSE, all.y=FALSE)

ibaqpy_ppb_aggregate$Tissues <-gsub("Esophageal epithelium", "Esophagus", ibaqpy_ppb_aggregate$Tissues)

DIA_DDA_ppb_merge <- merge(x=DDA_all_long, y=ibaqpy_ppb_aggregate,
                          by.x=c("Gene.ID","Gene.Symbol","Tissues"), by.y=c("EnsemblID", "Gene.Symbol","Tissues"),
                          al.x=FALSE, all.y=FALSE)
DIA_DDA_ppb_merge <- DIA_DDA_ppb_merge[DIA_DDA_ppb_merge$Gene.ID != "",]

ggplot(DIA_DDA_ppb_merge, aes(x=median_IbaqPpb, y=Median_ppb.iBAQ)) + 
  geom_point(size=0.5, alpha=0.2) + 
  scale_x_log10()+
  scale_y_log10()+
  xlab("Median ppb.iBAQ(DIA)")+
  ylab("Median ppb.iBAQ(DDA)")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x, size=0.4) +
  ggtitle("Protein abundance comparison\nppb.iBAQ(DDA-RefProt) v/s ibaqpy iBAQ(DIA-RefProt)")+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=12))+
  facet_wrap(~Tissues, scales = "free")

###################################################################
# 3. Compute iBAQ values from DIA-NN LFQ : divide by number of theoretical tryptic peptides
###################################################################
library("Biostrings")

if(Rootdata == "SwissProt"){
  fasta_file <- readAAStringSet("Human_OneProteinPerGeneSet_May2023_UP000005640_9606_PLUS_Contaminants.fasta")} else {
  fasta_file <- readAAStringSet("Human_UniProt_ReferenceProteome_withIsoforms_March2023_PLUS_Contaminants.fasta")
  }

fasta_seq <- as.data.frame(fasta_file)
fasta_seq <- tibble::rownames_to_column(fasta_seq, "name")
colnames(fasta_seq) <- c("name","seq")
fasta_seq$name <- gsub(" .*","", fasta_seq$name, perl=TRUE)
fasta_seq$number_of_tryptic_peptides <- "NA"

fragments <- list()

# Compute number of tryptic peptides in each sequence
# Look for Trypsin cleavage sites ----
for(i in 1:nrow(fasta_seq)){
  
  seq_id <-fasta_seq[i,1]
  seq <- fasta_seq[i,2]
  seq_length <- nchar(fasta_seq[i,2])
  
  print(seq_id)
  #get coordinates of Trypsin cleavage sites (Lysine or Arginine, but not followed by Proline)
  repl_seq <- seq
  repl_seq <- gsub("KP|RP","xx", repl_seq, perl=TRUE)
  seq_frag_coords <- data.frame(gregexpr(pattern ='K|R',repl_seq))
  fasta_seq[i, "number_of_tryptic_peptides"] <- nrow(seq_frag_coords)+1
}

fasta_seq$protein_id <-  gsub("sp\\||tr\\|", "", fasta_seq$name, perl=TRUE) 
fasta_seq$protein_id <-  gsub("\\|.*", "", fasta_seq$protein_id, perl=TRUE) 

#combined_ProteinLFQs_2ormore_pept_median_tissues


ProteinLFQs_2ormore_pept_Long_trypticpeptides <- merge(x=combined_ProteinLFQs_2ormore_pept_median_tissues, y=fasta_seq,
                                                  by.x=c("Protein.Group"), by.y=c("protein_id"),
                                                  all.x=FALSE, all.y=FALSE)

ProteinLFQs_2ormore_pept_Long_trypticpeptides <- ProteinLFQs_2ormore_pept_Long_trypticpeptides[complete.cases(ProteinLFQs_2ormore_pept_Long_trypticpeptides),]
ProteinLFQs_2ormore_pept_Long_trypticpeptides$Median_LFQ <- as.numeric(ProteinLFQs_2ormore_pept_Long_trypticpeptides$Median_LFQ)
ProteinLFQs_2ormore_pept_Long_trypticpeptides$number_of_tryptic_peptides <- as.numeric(ProteinLFQs_2ormore_pept_Long_trypticpeptides$number_of_tryptic_peptides)

ProteinLFQs_2ormore_pept_Long_trypticpeptides$ppb_LFQ_trypticpeptides <- ProteinLFQs_2ormore_pept_Long_trypticpeptides$Median_LFQ/ProteinLFQs_2ormore_pept_Long_trypticpeptides$number_of_tryptic_peptides


combined_DDA_DIA_ppb_trypticpeptides <- merge(x=ProteinLFQs_2ormore_pept_Long_trypticpeptides,
                                              y=DDA_all_long,
                                              by.x=c("Genes","Protein.Group","Tissues"), by.y=c("Gene.Symbol","Protein.Ids","Tissues"),
                                              all.x=FALSE, all.y=FALSE)

ggplot(combined_DDA_DIA_ppb_trypticpeptides, aes(x=ppb_LFQ_trypticpeptides, y=Median_ppb.iBAQ)) + 
  geom_point(size=0.5, alpha=0.2) + 
  scale_x_log10()+
  scale_y_log10()+
  xlab("Median iBAQ (DIA)")+
  ylab("Median ppb.iBAQ (DDA)")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x, size=0.4) +
  ggtitle("Protein abundance comparison\nppb.iBAQ(DDA-RefProt) v/s iBAQ.tryptic_peptides(DIA-RefProt 100KSet)")+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=12))+
  facet_wrap(~Tissues, scales = "free")
