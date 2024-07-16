#Script to add Gene Symbols to UniProt ID mapped table
#Andrew Collins' script generates mappings between UniPot_IDs and ENSGs, 
# but this does not have mappings to Gene Symbols.
# To add Gene Symbol information to this table:
#1, First map all UniProt_IDs from the table to Gene Symbols using UniProt ID mapper https://www.uniprot.org/id-mapping
# From database "UniProtKB AC/ID" To database "UniProtKB -> Gene Name"
#2. Then collapse multiple mappings to single row


uniprot_ids <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/idmapping_HUMAN.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

uniprot_genesymb <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/GeneSymbols_HUMAN.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

library(dplyr)

dups <- data.frame(table(uniprot_genesymb$UniProt_ID))
dups <- dups[dups$Freq > 1,]

foo <- uniprot_genesymb %>%
  group_by(UniProt_ID) %>%
  mutate(combined_symb = toString(GeneSymbol))

foo <- foo[-c(2)]
foo1 <- unique(foo)

colnames(foo1) <- c("UniProt_ID","GeneSymbol")
write.table(foo1, "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/GeneSymbols_collapsed_HUMAN.tsv", sep = "\t", row.names = FALSE, quote = FALSE )


########## PIG ###########
uniprot_ids_pig <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Chicken_Cow_Pig/ModelPig_idmapping.tsv", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

uniprot_pig_genesymb <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Chicken_Cow_Pig/GeneSymbols_PIG.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

dups <- data.frame(table(uniprot_pig_genesymb$UniProt_ID))
dups <- dups[dups$Freq > 1,]



model_pig_ids <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Chicken_Cow_Pig/ModelPig_idmapping.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
merged <- merge(x=uniprot_ids_pig, y=model_pig_ids,
                by.x=c("UniProt_Acc"), by.y=c("UniProt_Acc"),
                all.x=FALSE, all.y=FALSE)



