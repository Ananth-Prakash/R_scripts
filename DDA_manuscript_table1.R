tmp <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD015928/proteinGroups.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill =TRUE)

###### Get the numbers for publication
### Use this function for post-process
 postprocess <- function(tmp){

  tmp <- tmp[ tmp[, "Reverse"] != "+", ]
  tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]
  return(tmp)
 }

 tmp <- postprocess(tmp)
 
 protein_groups <- nrow(tmp)
 peptides <- sum(tmp$Peptides)
 unique_peptides <- sum(tmp$Unique.peptides)
 