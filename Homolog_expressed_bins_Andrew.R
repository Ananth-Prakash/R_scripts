# Heatmap showing homolog bin values in tissues across mouse human and rat

library(ggplot2)
library(reshape2)
library(viridis)
library(stringr)
library(dplyr)
library(plyr)
library(tidyr)
library(pagenum)
library(ggpubr)
library(gridExtra)

homolog_bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-bins-one_to_one_mappings_only_common_organs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# to check how many one-to-one mappings by counting number of times gene names have occurred in the dataset
# mapping count 1 in all three genes are one-to-one mappings and anything else are either one-to-many or many-to-many mappings
mappings_human <- ddply(homolog_bins,.(GeneName_human),nrow)
mappings_mouse <- ddply(homolog_bins,.(GeneName_mouse),nrow)
mappings_rat <- ddply(homolog_bins,.(GeneName_rat),nrow)

colnames(mappings_human) <- c("GeneName_human", "GeneName_counts_human")
colnames(mappings_mouse) <- c("GeneName_mouse", "GeneName_counts_mouse")
colnames(mappings_rat) <- c("GeneName_rat", "GeneName_counts_rat")

homolog_bins <- merge(homolog_bins, mappings_human,
                      by.x=c("GeneName_human"), by.y=c("GeneName_human"),
                      all.x=FALSE, all.y=FALSE)
homolog_bins <- merge(homolog_bins, mappings_mouse,
                      by.x=c("GeneName_mouse"), by.y=c("GeneName_mouse"),
                      all.x=FALSE, all.y=FALSE)
homolog_bins <- merge(homolog_bins, mappings_rat,
                      by.x=c("GeneName_rat"), by.y=c("GeneName_rat"),
                      all.x=FALSE, all.y=FALSE)

# Filter down to include only one-to-one mappings
homolog_bins <- homolog_bins[(homolog_bins$GeneName_counts_human ==1 & homolog_bins$GeneName_counts_mouse == 1 & homolog_bins$GeneName_counts_rat == 1),]

homolog_bins <- homolog_bins[order(homolog_bins$GeneName_human),]
homolog_bins$median <- apply(homolog_bins[,c(8:25)], 1, median, na.rm = T)

homolog_bins <- homolog_bins[homolog_bins$Present_in_samples != 0,]
data <- homolog_bins[,-c(1:7,26:30)]

data$Homologs <- factor(paste(homolog_bins$GeneName_human, homolog_bins$UniProt_ID_human, homolog_bins$mouse_ensg, homolog_bins$human_ensg, homolog_bins$rat_ensg, sep=" / "))


data_long <- gather(data, Organs, Bins, Brain_mouse:Testis_rat, factor_key=TRUE)
data_long$Species <- gsub(".*_","", data_long$Organs, perl=TRUE)
data_long$Organs <- gsub("_.*","", data_long$Organs, perl=TRUE)

data_long$Species <- gsub("mouse","Mouse",data_long$Species)
data_long$Species <- gsub("human","Human",data_long$Species)
data_long$Species <- gsub("rat","Rat",data_long$Species)


data_summary <- data
data_summary$Mouse_all_median <- apply(data_summary[,grepl("mouse", colnames(data_summary))],1,median, na.rm =TRUE)
data_summary$Human_all_median <- apply(data_summary[,grepl("human", colnames(data_summary))],1,median, na.rm =TRUE)
data_summary$Rat_all_median <- apply(data_summary[,grepl("rat", colnames(data_summary))],1,median, na.rm =TRUE)

data_summary$Mouse_count <- apply(data_summary[,grepl("_mouse", colnames(data_summary))],1,function(x) (length(which(!is.na(x)))))
data_summary$Human_count <- apply(data_summary[,grepl("_human", colnames(data_summary))],1,function(x) (length(which(!is.na(x)))))
data_summary$Rat_count <- apply(data_summary[,grepl("_rat", colnames(data_summary))],1,function(x) (length(which(!is.na(x)))))

# A subdata where orthologs are expressed is expressed in human, rat and mouse 
subdata <- data_summary[(data_summary$Human_count >= 3 & data_summary$Mouse_count >= 3  & data_summary$Rat_count >= 3 ),]

#take random 10 entries
subdata <- subdata[sample(nrow(subdata), 10), ]
subdata_long <- gather(subdata, Organs, Bins, Brain_human:Testis_rat, factor_key=TRUE)
subdata_long$Species <- gsub(".*_","", subdata_long$Organs, perl=TRUE)
subdata_long$Organs <- gsub("_.*","", subdata_long$Organs, perl=TRUE)

subdata_long$Species <- gsub("mouse","Mouse",subdata_long$Species)
subdata_long$Species <- gsub("human","Human",subdata_long$Species)
subdata_long$Species <- gsub("rat","Rat",subdata_long$Species)

subdata_long$Homologs <- gsub("\\/.*","", subdata_long$Homologs, perl=TRUE)
ggplot(subdata_long, aes(x=Organs,y=Homologs))+
  geom_tile(alpha=0, color = "gray",
            lwd = 0.3,
            linetype = 1)+
  geom_point(aes(colour = Species, 
                 size =Bins), position=position_dodge(width=0.5))+
  scale_y_discrete(limits=rev)+
  theme_bw()

x=0
#setPagenum(1)
input <- data.frame()
plot_list = list()

no_of_rows_per_plotpage <- 50
no_of_plot_pages <- ceiling(nrow(homolog_bins)/no_of_rows_per_plotpage)

for (i in 1:no_of_plot_pages) {
  
  subdata<-tail(head(data, n=50+x), n=50)
  
  input <- rbind(input,subdata)
  x=nrow(input)
  
  input_long <- gather(subdata, Organs, Bins, Brain_human:Testis_rat, factor_key=TRUE)
  input_long$Species <- gsub(".*_","", input_long$Organs, perl=TRUE)
  input_long$Organs <- gsub("_.*","", input_long$Organs, perl=TRUE)
  
  input_long$Species <- gsub("mouse","Mouse",input_long$Species)
  input_long$Species <- gsub("human","Human",input_long$Species)
  input_long$Species <- gsub("rat","Rat",input_long$Species)
  
  p = ggplot(input_long,aes(x=Organs,y=Homologs))+
    geom_tile(alpha=0, color = "gray",
              lwd = 0.3,
              linetype = 1)+
    geom_point(aes(colour = Species, 
                   size =Bins), position=position_dodge(width=0.5))+
    scale_y_discrete(limits=rev)+
    theme_bw()
  #pagenum(text="ABC Corp - ",x=.95, y=.2, just=c('right','bottom'))
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
pdf(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Plots/Homolog_binned_expression_comparison.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 15) # The height of the plot in inches
for (i in 1:no_of_plot_pages) {
  print(plot_list[[i]])
}

# Step 3: Run dev.off() to create the file!
dev.off()
