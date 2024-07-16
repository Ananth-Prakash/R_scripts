library(dplyr)

Species <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/Dataset_lists_Species.txt', quote = "\"", sep="\t", header=TRUE)
PubMed <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/Dataset_lists_PubMed.txt', quote = "\"", sep="\t", header=TRUE)
Size <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/Dataset_lists_FileSize.txt', quote = "\"", sep="\t", header=TRUE)
DIA_accessions <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/swath_potentially_reusable_Accessions.txt', quote = "\"", sep="\t", header=TRUE)



Merge1 <- merge(x=Species, y= PubMed,
                by.x=c("ACCESSION"), by.y=c("ACCESSION"),
                all.x=TRUE, all.y=TRUE)

Merge2 <- merge(x=Merge1, y= Size,
                by.x=c("ACCESSION"), by.y=c("ACCESSION"),
                all.x=TRUE, all.y=TRUE)

Merge3 <- merge(x=Merge2, y=DIA_accessions,
                by.x=c("ACCESSION"), by.y=c("accession"),
                all.x=FALSE, all.y=FALSE)


Aggregated_df <- Merge3 %>%
  group_by(ACCESSION, TITLE, PROJECT_DESCRIPTION, SAMPLE_PROC, DATA_PROC, SPECIES, EMAIL,
           FIRST_NAME, LAST_NAME, FILE_SIZE, COUNT) %>%
  summarise(PUBMED_ID = toString(PUBMED_ID))

Human_samples <- Aggregated_df[Aggregated_df$SPECIES == "Homo sapiens (Human)",]

write.table(Human_samples, "/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/SWATH-DIA_datasets_Human_annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE )
write.table(Aggregated_df, "/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/SWATH-DIA_datasets_All_species_annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE )




annotated <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/SWATH-DIA_datasets_All_species_annotated.txt', quote = "\"", sep="\t", header=TRUE)
instrument <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/instrument-DAVID1.txt', quote = "\"", sep="\t", header=TRUE)
tissue <- read.table('/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/tissue-DAVID1.txt', quote = "\"", sep="\t", header=TRUE)


merge4 <- merge(x=annotated, y=tissue,
                by.x=c("ACCESSION"), by.y=c("ACCESSION"),
                all.x=TRUE, all.y=FALSE)
merge5 <- merge(x=merge4, y=instrument,
               by.x=c("ACCESSION"), by.y=c("ACCESSION"),
               all.x=TRUE, all.y=FALSE)

Aggregated_1_df <- merge5 %>%
  group_by(ACCESSION, TITLE, PROJECT_DESCRIPTION, SAMPLE_PROC, DATA_PROC, SPECIES, EMAIL,
           FIRST_NAME, LAST_NAME, FILE_SIZE, COUNT) %>%
  summarise(PUBMED_ID = toString(PUBMED_ID),
            NAME = toString(NAME),
            INSTRUMENT = toString(INSTRUMENT))

write.table(Aggregated_1_df, "/Users/ananth/Documents/MaxQuant_Bechmarking/SWATH-DIA/SWATH-DIA_datasets_All_species_tissue_instrument_annotated.txt", sep = "\t", row.names = FALSE, quote = FALSE )


