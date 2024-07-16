library(dplyr)
library(stringr)
# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")

library(mygene)

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/")

Human_file_list <- data.frame(Dataset = c(
# "PXD010154_Ananth", 
# "PXD005819_33threads_yoda", 
# "PXD004143",
# "PXD006233",
# "PXD012755", 
# "PXD001608_30threads_yoda", 
# "PXD002029", 
# "PXD000547", 
# "PXD000548", 
# "PXD010271", 
# "PXD004332",
# "PXD006675",
# "PXD008934", 
# "Synapse-AD/ACT_DorsoLateralPreFrontalCortex", 
# "Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex",
# "Synapse-AD/Banner_DorsoLateralPreFrontalCortex",
# "Synapse-AD/BLSA_DorsoLateralPreFrontalCortex",
# "Synapse-AD/BLSA_Precuneus", 
# "Synapse-AD/Mayo_TemporalCortex", 
# "Synapse-AD/MountSinai_FrontalPole",
# "Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases",
# "PXD012131",
# "PXD020187", 
# "PXD015079"
))

Rat_file_list <- data.frame(Dataset = c(
#  "PXD012677",
#  "PXD006692",
#  "PXD016793",
#  "PXD004364",
#  "PXD001839",
#  "PXD013543",
#  "PXD016958",
#  "PXD003375",
#  "PXD015928"
))

Mouse_file_list <- data.frame(Dataset = c(
# "PXD000867",
# "PXD000288",
# "PXD003155",
# "PXD004612",
# "PXD005230",
# "PXD009909",
# "PXD012307",
# "PXD009639",
# "PXD019394",
# "PXD012636",
# "PXD019431",
# "PXD022614",
# "PXD004496",
  "PXD008736"
  ))


proteinGroups_decoy_list <- list()

file_list <- Mouse_file_list

for(i in 1: nrow(file_list)){
  
  datasetID <- file_list[i,]
  print(datasetID)
#  tmp  <- read.table(file = paste("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/", datasetID, "/proteinGroups.txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
#  tmp  <- read.table(file=paste("/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/", datasetID, "/proteinGroups.txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  tmp  <- read.table(file=paste("/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/", datasetID, "/proteinGroups.txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  
    print(nrow(tmp))
    proteinGroups_decoy_list[[i]] <- tmp[ tmp[, "Reverse"] == "+", c("Majority.protein.IDs"), drop=FALSE]
}

all_decoys <- bind_rows(proteinGroups_decoy_list)

all_decoys <- unique(all_decoys)

all_decoys$"Decoy.ENSG" <- "NA"
all_decoys$"Decoy.Gene.Symbol" <- "NA"
all_decoys$"unique.gene.count" <- "NA"


for(i in 1:nrow(all_decoys)){
  
  x <- data.frame(strsplit(all_decoys[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
  x_temp <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  x_temp <- gsub("REV__", "", x_temp, perl=TRUE)
  #x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  #x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
  f = file()
  sink(file=f)
  #res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("human")),
  #                error = function(e) {print(0)})
  #res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("rat")),
  #                error = function(e) {print(0)})
  res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("mouse")),
                  error = function(e) {print(0)})
  
  # Covid2 taxonomy id: 2697049
  sink()
  close(f)
  
  if (class(res)=="DFrame" | class(res) == "DataFrame"){
    all_decoys[ i, "Decoy.ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    #if(all_decoys[i,"ENSG"] == ""){
    #  all_decoys[ i, "ENSG"] <- paste( unique(unlist(res$ensembl[!is.na(res$ensembl)])), collapse = ";")
    #}
    all_decoys[ i, "Decoy.Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")}
  temp_symb <- all_decoys[i,"Decoy.Gene.Symbol"]
  all_decoys[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
  
  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(all_decoys),1)),"%"))
}

all_decoys <- all_decoys[all_decoys$Decoy.ENSG != "" , ]
all_decoys <- all_decoys[all_decoys$Decoy.ENSG != "NA" , ]
all_decoys <- all_decoys[all_decoys$Decoy.Gene.Symbol != "NA" , ]
all_decoys <- all_decoys[all_decoys$Decoy.Gene.Symbol != "" , ]
nrow(unique(data.frame(all_decoys$Decoy.ENSG)))
all_decoys <- all_decoys[all_decoys$unique.gene.count == 1,] ##### comment this later
nrow(unique(data.frame(all_decoys$Decoy.ENSG)))

all_decoys_mapped <- all_decoys

get_decoy <- function(datasetID){
   tmp  <- read.table(file = paste("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/", datasetID, "/proteinGroups-tissue_names.txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
   decoy <- tmp[ tmp[, "Reverse"] == "+", c("Majority.protein.IDs"), drop=FALSE]

   decoy_merged <- merge(x=all_decoys_mapped, y=decoy,
                         by.x=c("Majority.protein.IDs"), by.y=c("Majority.protein.IDs"),
                         all.x=FALSE, all.y=FALSE)

   decoy_merged <- decoy_merged[,c("Decoy.Gene.Symbol"), drop=FALSE]
   decoy_merged$x <- decoy_merged$Decoy.Gene.Symbol
   colnames(decoy_merged) <- c("Decoy.Gene.Symbol",paste("Decoy.Gene.Symbol",datasetID,sep="_"))
   return(decoy_merged)
}

merge_decoys <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Decoy.Gene.Symbol"), by.y=c("Decoy.Gene.Symbol"), all.x=TRUE, all.y=TRUE)
  return(merged)
}

data1_decoy <- get_decoy("PXD010154_Ananth")
data2_decoy <- get_decoy("PXD005819_33threads_yoda")
data3_decoy <- get_decoy("PXD004143")
data4_decoy <- get_decoy("PXD006233")
data5_decoy <- get_decoy("PXD012755")
data6_decoy <- get_decoy("PXD001608_30threads_yoda")
data7_decoy <- get_decoy("PXD002029")
data8_decoy <- get_decoy("PXD000547")
data9_decoy <- get_decoy("PXD000548")
data10_decoy <- get_decoy("PXD010271")
data11_decoy <- get_decoy("PXD004332")
data12_decoy <- get_decoy("PXD006675")
data13_decoy <- get_decoy("PXD008934")
data14_decoy <- get_decoy("Synapse-AD/ACT_DorsoLateralPreFrontalCortex")
data15_decoy <- get_decoy("Synapse-AD/Aging_Baltimore-JohnsHopkins_DorsoLateralPreFrontalCortex")
data16_decoy <- get_decoy("Synapse-AD/Banner_DorsoLateralPreFrontalCortex")
data17_decoy <- get_decoy("Synapse-AD/BLSA_DorsoLateralPreFrontalCortex")
data18_decoy <- get_decoy("Synapse-AD/BLSA_Precuneus")
data19_decoy <- get_decoy("Synapse-AD/Mayo_TemporalCortex")
data20_decoy <- get_decoy("Synapse-AD/MountSinai_FrontalPole")
data21_decoy <- get_decoy("Synapse-AD/UniPenn_MultipleNeurodegenerativeDiseases")
data22_decoy <- get_decoy("PXD012131")
data23_decoy <- get_decoy("PXD020187")
data24_decoy <- get_decoy("PXD015079")

all_data_decoys_merged <- merge_decoys(all_decoys_mapped, data1_decoy)
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data2_decoy)
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data3_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data4_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data5_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data6_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data7_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data8_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data9_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data10_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data11_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data12_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data13_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data14_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data15_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data16_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data17_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data18_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data19_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data20_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data21_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data22_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data23_decoy) 
all_data_decoys_merged <- merge_decoys(all_data_decoys_merged, data24_decoy) 


all_data_decoys_merged <- cbind(all_data_decoys_merged, data.frame(Present_in_number_of_samples = apply(all_data_decoys_merged[5:ncol(all_data_decoys_merged)], 1, function(x) length(which(x != "NaN")) )))
all_data_decoys_merged <- all_data_decoys_merged[,-c(2)]
all_data_decoys_merged = all_data_decoys_merged[!duplicated(all_data_decoys_merged$Decoy.Gene.Symbol),]

write.table(all_data_decoys_merged, file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Decoy_as_individuals_distribution_in_all_datasets.txt", sep = "\t", row.names = FALSE, quote = FALSE )


ggplot(all_data_decoys_merged, aes(x=Present_in_number_of_samples)) + 
  geom_histogram( aes(fill="red") )+ 
  xlab("Spread across number of datasets")+
  ylab("Reverse decoys")+
  theme_bw()+
  theme(legend.position = "none")+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.7) +
  ggtitle("Distribution of reverse decoys (n=1735) collected from all datasets (n=24)\n mapped to more than one canonical gene")



  