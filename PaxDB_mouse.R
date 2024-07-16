# Compare mouse PaxDB data with MaxQuant
library(mygene)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

PaxDB_brain  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_BRAIN_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_brain) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_brain$Ensembl_Protein <- gsub(".*\\.","", PaxDB_brain$Ensembl_Protein, perl=TRUE)

PaxDB_adipose  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_BROWN_ADIPOSE_TISSUE_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_adipose) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_adipose$Ensembl_Protein <- gsub(".*\\.","", PaxDB_adipose$Ensembl_Protein, perl=TRUE)

PaxDB_heart  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_HEART_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_heart) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_heart$Ensembl_Protein <- gsub(".*\\.","", PaxDB_heart$Ensembl_Protein, perl=TRUE)

PaxDB_kidney  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_KIDNEY_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_kidney) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_kidney$Ensembl_Protein <- gsub(".*\\.","", PaxDB_kidney$Ensembl_Protein, perl=TRUE)

PaxDB_liver  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_LIVER_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_liver) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_liver$Ensembl_Protein <- gsub(".*\\.","", PaxDB_liver$Ensembl_Protein, perl=TRUE)

PaxDB_lung  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_LUNG_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_lung) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_lung$Ensembl_Protein <- gsub(".*\\.","", PaxDB_lung$Ensembl_Protein, perl=TRUE)

PaxDB_pancreas  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_PANCREAS_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_pancreas) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_pancreas$Ensembl_Protein <- gsub(".*\\.","", PaxDB_pancreas$Ensembl_Protein, perl=TRUE)

PaxDB_spleen  <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/PaxDB/PaxDB_Mouse_SPLEEN_10090-integrated.txt" , quote = "\"", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(PaxDB_spleen) <- c("InternalID","Ensembl_Protein","PaxDB_abundance")
PaxDB_spleen$Ensembl_Protein <- gsub(".*\\.","", PaxDB_spleen$Ensembl_Protein, perl=TRUE)


##### MaxQuant analysis

MaxQuant_mouse <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/Gene_distribution_in_organs-GeneNames-Median_intensities-Mouse.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
MaxQuant_mouse <- MaxQuant_mouse[,-c(15:18)]


data.to.map <- MaxQuant_mouse
data.to.map$"ENSP" <- "NA"
data.to.map$"ENSG" <- "NA"


for(i in 1:nrow(data.to.map)){
  
  x <- data.to.map[ i, "GeneName"]
  print(x)
  #x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
  #x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  #x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
  f = file()
  sink(file=f)
  res <- tryCatch(queryMany(x, scopes="symbol", fields=c("ensembl.gene","ensembl.protein", "symbol"), species=c("mouse")),
                error = function(e) {print(0)})
  # Covid2 taxonomy id: 2697049
  sink()
  close(f)

  if (class(res)=="DFrame" | class(res) == "DataFrame"){
    data.to.map[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    data.to.map[ i, "ENSP"] <- paste( unique(unlist(res$ensembl.protein[!is.na(res$ensembl.protein)])), collapse = ";")
    #if(data.to.map[i,"ENSG"] == ""){
    #  data.to.map[ i, "ENSG"] <- paste( unique(unlist(res$ensembl[!is.na(res$ensembl)])), collapse = ";")
    #}
     #temp_symb <- data.to.map[i,"ENSP"]
  #data.to.map[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1

  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(data.to.map),1)),"%"))
  }
}

backup <- data.to.map

#splitting ENSP into individual entries  
data.to.map <- data.to.map %>% 
  mutate(ENSP = strsplit(as.character(ENSP), ";")) %>% 
  unnest(ENSP)
data.to.map <- unique(data.to.map)

MaxQuant_brain <- data.to.map[,c("ENSP","Brain_median")]
Merged_brain <- merge(x=PaxDB_brain, y=MaxQuant_brain,
                      by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                      all.x=FALSE, all.y=FALSE)
Merged_brain$PaxDB_abundance_log2 <- log2(Merged_brain$PaxDB_abundance)
Merged_brain$MaxQuant_abundance_log2 <- log2(Merged_brain$Brain_median)
Merged_brain$Organ <- rep("Brain", nrow(Merged_brain))

MaxQuant_brain <- data.to.map[,c("ENSP","Brain_median")]


MaxQuant_heart <- data.to.map[,c("ENSP","Heart_median")]
Merged_heart <- merge(x=PaxDB_heart, y=MaxQuant_heart,
                      by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                      all.x=FALSE, all.y=FALSE)
Merged_heart$PaxDB_abundance_log2 <- log2(Merged_heart$PaxDB_abundance)
Merged_heart$MaxQuant_abundance_log2 <- log2(Merged_heart$Heart_median)
Merged_heart$Organ <- rep("Heart", nrow(Merged_heart))

MaxQuant_kidney <- data.to.map[,c("ENSP","Kidney_median")]
Merged_kidney <- merge(x=PaxDB_kidney, y=MaxQuant_kidney,
                      by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                      all.x=FALSE, all.y=FALSE)
Merged_kidney$PaxDB_abundance_log2 <- log2(Merged_kidney$PaxDB_abundance)
Merged_kidney$MaxQuant_abundance_log2 <- log2(Merged_kidney$Kidney_median)
Merged_kidney$Organ <- rep("Kidney", nrow(Merged_kidney))

MaxQuant_liver <- data.to.map[,c("ENSP","Liver_median")]
Merged_liver <- merge(x=PaxDB_liver, y=MaxQuant_liver,
                       by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                       all.x=FALSE, all.y=FALSE)
Merged_liver$PaxDB_abundance_log2 <- log2(Merged_liver$PaxDB_abundance)
Merged_liver$MaxQuant_abundance_log2 <- log2(Merged_liver$Liver_median)
Merged_liver$Organ <- rep("Liver", nrow(Merged_liver))

MaxQuant_lung <- data.to.map[,c("ENSP","Lung_median")]
Merged_lung <- merge(x=PaxDB_lung, y=MaxQuant_lung,
                      by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                      all.x=FALSE, all.y=FALSE)
Merged_lung$PaxDB_abundance_log2 <- log2(Merged_lung$PaxDB_abundance)
Merged_lung$MaxQuant_abundance_log2 <- log2(Merged_lung$Lung_median)
Merged_lung$Organ <- rep("Lung", nrow(Merged_lung))

MaxQuant_pancreas <- data.to.map[,c("ENSP","Pancreas_median")]
Merged_pancreas <- merge(x=PaxDB_pancreas, y=MaxQuant_pancreas,
                     by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                     all.x=FALSE, all.y=FALSE)
Merged_pancreas$PaxDB_abundance_log2 <- log2(Merged_pancreas$PaxDB_abundance)
Merged_pancreas$MaxQuant_abundance_log2 <- log2(Merged_pancreas$Pancreas_median)
Merged_pancreas$Organ <- rep("Pancreas", nrow(Merged_pancreas))

MaxQuant_spleen <- data.to.map[,c("ENSP","Spleen_median")]
Merged_spleen <- merge(x=PaxDB_spleen, y=MaxQuant_spleen,
                         by.x=c("Ensembl_Protein"), by.y=c("ENSP"),
                         all.x=FALSE, all.y=FALSE)
Merged_spleen$PaxDB_abundance_log2 <- log2(Merged_spleen$PaxDB_abundance)
Merged_spleen$MaxQuant_abundance_log2 <- log2(Merged_spleen$Spleen_median)
Merged_spleen$Organ <- rep("Spleen", nrow(Merged_spleen))


plotdata <- rbind(Merged_brain[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")],
                 Merged_heart[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")],
                 Merged_kidney[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")],
                 Merged_liver[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")],
                 Merged_lung[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")],
                 Merged_pancreas[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")],
                 Merged_spleen[,c("Ensembl_Protein","PaxDB_abundance_log2","MaxQuant_abundance_log2","Organ")])

plotdata$PaxDB_abundance_log2[plotdata$PaxDB_abundance_log2 == "-Inf"] <- NA
plotdata$MaxQuant_abundance_log2[plotdata$MaxQuant_abundance_log2 == "-Inf"] <- NA

ggplot(plotdata, aes(x=MaxQuant_abundance_log2, y=PaxDB_abundance_log2)) + 
  geom_point(size=0.05, alpha=0.5) + 
  xlab("log2(FOT normalised iBAQ)  Mouse_Wang et al 2021")+
  ylab("log2(abundance)  Mouse_PaxDB")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse protein abundance comparison (Wang et al 2021 vs PaxDB)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  facet_wrap(~Organ)




