foo1 <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/PXD004143/proteinGroups.txt" , quote = "\"", check.names=FALSE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

get_subdata <- function(tmp){
  tmp <- tmp[ tmp[, "Reverse"] != "+", ]
  tmp <- tmp[ tmp[, "Potential contaminant"] != "+", ]
  ## consider those protein groups with more than 1 peptides matched
  #tmp <- tmp[ tmp[, "Peptides"] > 1, ]
  
  #tmp <- tmp[,c(2,10,12,127,170,171, 326, 248,288)] # For PXD010154
  tmp <- tmp[,c("Majority protein IDs","Peptides", "Unique peptides","Sequence coverage [%]","Q-value","Score","MS/MS count","Intensity","iBAQ")]
  #tmp <- tmp[,c("Majority.protein.IDs","Peptides", "Unique.peptides","Sequence.coverage....","Q.value","Score","MS.MS.count","Intensity","iBAQ")]
   return(tmp)
}

subfoo1 <- get_subdata(foo1)

write.table(subfoo1, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/PXD010154_subdata.txt", sep = "\t", row.names = FALSE, quote = FALSE )
