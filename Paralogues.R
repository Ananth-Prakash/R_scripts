#https://support.bioconductor.org/p/105070/

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
rat   = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/')

#####################
#1. Human
#####################
#gene_id <- "ENSG00000165156"
human_gene_id <- read.table(file = "human_genes_for_paralogs_search.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
human_gene_id <- data.frame(human_gene_id[!duplicated(human_gene_id), ])

human_paralogs_list <- list()
#If the loop crashes, restart from the last iteration number in the list
for(i in 1:nrow(human_gene_id)){
  human_paralogs_list[[i]] <- getBM(attributes = c("ensembl_gene_id", 
                                "external_gene_name",
                                "hsapiens_paralog_ensembl_gene", 
                                "hsapiens_paralog_associated_gene_name"),
                 #for "ENSG00000165156"
                 filters = "ensembl_gene_id", 
                 #for "ZHX1"
                 #filters = "external_gene_name", 
                 values = human_gene_id[i,1],
                 mart = human)
  print(paste0("Finding paralogues for ", human_gene_id[i,1], "... ", as.character(round(i*100/nrow(human_gene_id),1)),"%"))
}
backup_human <- human_paralogs_list
human_paralogs_list <- Filter(function(x) nrow(x) > 0, human_paralogs_list)

human_paralogs <- bind_rows(human_paralogs_list)

human_paralogs <- human_paralogs[complete.cases(human_paralogs), ]
human_paralogs <- human_paralogs[human_paralogs$hsapiens_paralog_associated_gene_name != "",]

colnames(human_paralogs) <- c("gene_id","gene_name","paralog_id","paralog_name")
foo <- human_paralogs

#remove duplicated symmetrical mapping 
#(ie., remove one instance of a duplicated row where gene-paralog is the same as paralog-gene)
#(ex. ENSG00000103653	CSK ENSG00000000938	FGR is the same mapping as
#     ENSG00000000938	FGR	ENSG00000103653	CSK, therefore keep only one instance of this symmetrical mapping)

for(i in 1:nrow(human_paralogs)){
  gene <- human_paralogs[i,"gene_name"]
  paralog <- human_paralogs[i,"paralog_name"]
  if(!is.na(gene) | !is.na(paralog)){
    human_paralogs <- subset(human_paralogs, gene_name != paralog | paralog_name != gene)
  }
}

write.table(human_paralogs, file = "Backup_mapped_paralogues-human.txt", sep = "\t", row.names = FALSE, quote = FALSE )

#homolog_ppb <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-ppb-one_to_one_mappings_only_common_organs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
homolog_ppb <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-ppb.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

human_paralogs_left <- merge(x=human_paralogs, y=homolog_ppb[grepl("human_ensg|median_human", colnames(homolog_ppb), ignore.case = TRUE)],
                        by.x=c("gene_id"), by.y=c("human_ensg"),
                        all.x=FALSE, all.y=FALSE)

human_paralogs_right <- merge(x=human_paralogs, y=homolog_ppb[grepl("human_ensg|median_human", colnames(homolog_ppb), ignore.case = TRUE)],
                             by.x=c("paralog_id"), by.y=c("human_ensg"),
                             all.x=FALSE, all.y=FALSE)

human_left <- human_paralogs_left[grepl("gene_|_median_human",colnames(human_paralogs_left), ignore.case = TRUE)]
human_right <- human_paralogs_right[grepl("paralog_|_median_human",colnames(human_paralogs_right), ignore.case = TRUE)]

colnames(human_left) <- gsub("_median_human","_gene",colnames(human_left), perl=TRUE)
colnames(human_right) <- gsub("_median_human","_paralog",colnames(human_right), perl=TRUE)

human_paralogs_merged <- merge(x=human_paralogs, y=human_left,
                              by.x=c("gene_id","gene_name"),
                              by.y=c("gene_id","gene_name"), all.x=FALSE, all.y=FALSE)

human_paralogs_merged <- merge(x=human_paralogs_merged, y=human_right,
                              by.x=c("paralog_id","paralog_name"),
                              by.y=c("paralog_id","paralog_name"), all.x=FALSE, all.y=FALSE)

human_paralogs_merged <- human_paralogs_merged[!duplicated(human_paralogs_merged), ]

write.table(human_paralogs_merged, file = "Paralogues_mapping_ppb_expression_organs-human.txt", sep = "\t", row.names = FALSE, quote = FALSE )

human_long <- gather(human_paralogs_merged, condition, expression, AdiposeTissue_gene:VermiformAppendix_paralog, factor_key=TRUE)
human_long$Organ <- gsub("_.*","",human_long$condition, perl=TRUE)
human_long$Type <- gsub(".*_","",human_long$condition, perl=TRUE)
human_long$Species <- rep("Human",nrow(human_long))
human_long <- subset(human_long, select=-c(condition))

human_long_gene <- human_long[human_long$Type == "gene",]
human_long_paralog <- human_long[human_long$Type == "paralog",]

colnames(human_long_gene)[which(names(human_long_gene) == "expression")] <- "Gene"
colnames(human_long_paralog)[which(names(human_long_paralog) == "expression")] <- "Paralog"

human_long_all <- merge(x=human_long_gene, y=human_long_paralog,
                        by.x=c("paralog_id","paralog_name","gene_id","gene_name","Organ","Species"),
                        by.y=c("paralog_id","paralog_name","gene_id","gene_name","Organ","Species"),
                        all.x=FALSE, all.y=FALSE)
human_long_all <- subset(human_long_all, select=-c(Type.x, Type.y))
human_long_all <- human_long_all[order(human_long_all$Organ),]

write.table(human_long_all, file = "Paralogues_mapping_ppb_expression_organs-human-longformat.txt", sep = "\t", row.names = FALSE, quote = FALSE )


#####################
#2. Mouse
#####################
mouse_gene_id <- read.table(file = "mouse_genes_for_paralogs_search.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_gene_id <- data.frame(mouse_gene_id[!duplicated(mouse_gene_id), ])

mouse_paralogs_list <- list()
#If the loop crashes, restart from the last iteration number in the list
for(i in 1:nrow(mouse_gene_id)){
  mouse_paralogs_list[[i]] <- getBM(attributes = c("ensembl_gene_id", 
                                                   "external_gene_name",
                                                   "mmusculus_paralog_ensembl_gene", 
                                                   "mmusculus_paralog_associated_gene_name"),
                                    #for "ENSG00000165156"
                                    filters = "ensembl_gene_id", 
                                    #for "ZHX1"
                                    #filters = "external_gene_name", 
                                    values = mouse_gene_id[i,1],
                                    mart = mouse)
  print(paste0("Finding paralogues for ", mouse_gene_id[i,1], "... ", as.character(round(i*100/nrow(mouse_gene_id),1)),"%"))
}
backup_mouse <- mouse_paralogs_list
#Remove rows values in list that have 0 rows
mouse_paralogs_list <- Filter(function(x) nrow(x) > 0, mouse_paralogs_list)

mouse_paralogs <- bind_rows(mouse_paralogs_list)

mouse_paralogs <- mouse_paralogs[complete.cases(mouse_paralogs), ]
mouse_paralogs <- mouse_paralogs[mouse_paralogs$mmusculus_paralog_associated_gene_name != "",]

colnames(mouse_paralogs) <- c("gene_id","gene_name","paralog_id","paralog_name")
foo <- mouse_paralogs

#remove duplicated symmetrical mapping 
#(ie., remove one instance of a duplicated row where gene-paralog is the same as paralog-gene)
#(ex. ENSG00000103653	CSK ENSG00000000938	FGR is the same mapping as
#     ENSG00000000938	FGR	ENSG00000103653	CSK, therefore keep only one instance of this symmetrical mapping)

for(i in 1:nrow(mouse_paralogs)){
  gene <- mouse_paralogs[i,"gene_name"]
  paralog <- mouse_paralogs[i,"paralog_name"]
  if(!is.na(gene) | !is.na(paralog)){
    mouse_paralogs <- subset(mouse_paralogs, gene_name != paralog | paralog_name != gene)
  }
}

write.table(mouse_paralogs, file = "Backup_mapped_paralogues-mouse.txt", sep = "\t", row.names = FALSE, quote = FALSE )

homolog_ppb <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-ppb.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

mouse_paralogs_left <- merge(x=mouse_paralogs, y=homolog_ppb[grepl("mouse_ensg|median_mouse", colnames(homolog_ppb), ignore.case = TRUE)],
                             by.x=c("gene_id"), by.y=c("mouse_ensg"),
                             all.x=FALSE, all.y=FALSE)

mouse_paralogs_right <- merge(x=mouse_paralogs, y=homolog_ppb[grepl("mouse_ensg|median_mouse", colnames(homolog_ppb), ignore.case = TRUE)],
                              by.x=c("paralog_id"), by.y=c("mouse_ensg"),
                              all.x=FALSE, all.y=FALSE)

mouse_left <- mouse_paralogs_left[grepl("gene_|_median_mouse",colnames(mouse_paralogs_left), ignore.case = TRUE)]
mouse_right <- mouse_paralogs_right[grepl("paralog_|_median_mouse",colnames(mouse_paralogs_right), ignore.case = TRUE)]

colnames(mouse_left) <- gsub("_median_mouse","_gene",colnames(mouse_left), perl=TRUE)
colnames(mouse_right) <- gsub("_median_mouse","_paralog",colnames(mouse_right), perl=TRUE)

mouse_paralogs_merged <- merge(x=mouse_paralogs, y=mouse_left,
                               by.x=c("gene_id","gene_name"),
                               by.y=c("gene_id","gene_name"), all.x=FALSE, all.y=FALSE)

mouse_paralogs_merged <- merge(x=mouse_paralogs_merged, y=mouse_right,
                               by.x=c("paralog_id","paralog_name"),
                               by.y=c("paralog_id","paralog_name"), all.x=FALSE, all.y=FALSE)

mouse_paralogs_merged <- mouse_paralogs_merged[!duplicated(mouse_paralogs_merged), ]

write.table(mouse_paralogs_merged, file = "Paralogues_mapping_ppb_expression_organs-mouse.txt", sep = "\t", row.names = FALSE, quote = FALSE )

mouse_long <- gather(mouse_paralogs_merged, condition, expression, Brain_gene:ArticularCartilage_paralog, factor_key=TRUE)
mouse_long$Organ <- gsub("_.*","",mouse_long$condition, perl=TRUE)
mouse_long$Type <- gsub(".*_","",mouse_long$condition, perl=TRUE)
mouse_long$Species <- rep("Mouse",nrow(mouse_long))
mouse_long <- subset(mouse_long, select=-c(condition))

mouse_long_gene <- mouse_long[mouse_long$Type == "gene",]
mouse_long_paralog <- mouse_long[mouse_long$Type == "paralog",]

colnames(mouse_long_gene)[which(names(mouse_long_gene) == "expression")] <- "Gene"
colnames(mouse_long_paralog)[which(names(mouse_long_paralog) == "expression")] <- "Paralog"

mouse_long_all <- merge(x=mouse_long_gene, y=mouse_long_paralog,
                        by.x=c("paralog_id","paralog_name","gene_id","gene_name","Organ","Species"),
                        by.y=c("paralog_id","paralog_name","gene_id","gene_name","Organ","Species"),
                        all.x=FALSE, all.y=FALSE)
mouse_long_all <- subset(mouse_long_all, select=-c(Type.x, Type.y))
mouse_long_all <- mouse_long_all[order(mouse_long_all$Organ),]

write.table(mouse_long_all, file = "Paralogues_mapping_ppb_expression_organs-mouse-longformat.txt", sep = "\t", row.names = FALSE, quote = FALSE )


#####################
#3. Rat
#####################
rat_gene_id <- read.table(file = "rat_genes_for_paralogs_search.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
rat_gene_id <- data.frame(rat_gene_id[!duplicated(rat_gene_id), ])

rat_paralogs_list <- list()
#If the loop crashes, restart from the last iteration number in the list
for(i in 1:nrow(rat_gene_id)){
  rat_paralogs_list[[i]] <- getBM(attributes = c("ensembl_gene_id", 
                                                   "external_gene_name",
                                                   "rnorvegicus_paralog_ensembl_gene", 
                                                   "rnorvegicus_paralog_associated_gene_name"),
                                    #for "ENSG00000165156"
                                    filters = "ensembl_gene_id", 
                                    #for "ZHX1"
                                    #filters = "external_gene_name", 
                                    values = rat_gene_id[i,1],
                                    mart = rat)
  print(paste0("Finding paralogues for ", rat_gene_id[i,1], "... ", as.character(round(i*100/nrow(rat_gene_id),1)),"%"))
}
backup_rat <- rat_paralogs_list
#Remove rows values in list that have 0 rows
rat_paralogs_list <- Filter(function(x) nrow(x) > 0, rat_paralogs_list)

rat_paralogs <- bind_rows(rat_paralogs_list)

rat_paralogs <- rat_paralogs[complete.cases(rat_paralogs), ]
rat_paralogs <- rat_paralogs[rat_paralogs$rnorvegicus_paralog_associated_gene_name != "",]

colnames(rat_paralogs) <- c("gene_id","gene_name","paralog_id","paralog_name")
foo <- rat_paralogs

#remove duplicated symmetrical mapping 
#(ie., remove one instance of a duplicated row where gene-paralog is the same as paralog-gene)
#(ex. ENSG00000103653	CSK ENSG00000000938	FGR is the same mapping as
#     ENSG00000000938	FGR	ENSG00000103653	CSK, therefore keep only one instance of this symmetrical mapping)

for(i in 1:nrow(rat_paralogs)){
  gene <- rat_paralogs[i,"gene_name"]
  paralog <- rat_paralogs[i,"paralog_name"]
  if(!is.na(gene) | !is.na(paralog)){
    rat_paralogs <- subset(rat_paralogs, gene_name != paralog | paralog_name != gene)
  }
}

write.table(rat_paralogs, file = "Backup_mapped_paralogues-rat.txt", sep = "\t", row.names = FALSE, quote = FALSE )

homolog_ppb <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-ppb.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

rat_paralogs_left <- merge(x=rat_paralogs, y=homolog_ppb[grepl("rat_ensg|median_rat", colnames(homolog_ppb), ignore.case = TRUE)],
                             by.x=c("gene_id"), by.y=c("rat_ensg"),
                             all.x=FALSE, all.y=FALSE)

rat_paralogs_right <- merge(x=rat_paralogs, y=homolog_ppb[grepl("rat_ensg|median_rat", colnames(homolog_ppb), ignore.case = TRUE)],
                              by.x=c("paralog_id"), by.y=c("rat_ensg"),
                              all.x=FALSE, all.y=FALSE)

rat_left <- rat_paralogs_left[grepl("gene_|_median_rat",colnames(rat_paralogs_left), ignore.case = TRUE)]
rat_right <- rat_paralogs_right[grepl("paralog_|_median_rat",colnames(rat_paralogs_right), ignore.case = TRUE)]

colnames(rat_left) <- gsub("_median_rat","_gene",colnames(rat_left), perl=TRUE)
colnames(rat_right) <- gsub("_median_rat","_paralog",colnames(rat_right), perl=TRUE)

rat_paralogs_merged <- merge(x=rat_paralogs, y=rat_left,
                               by.x=c("gene_id","gene_name"),
                               by.y=c("gene_id","gene_name"), all.x=FALSE, all.y=FALSE)

rat_paralogs_merged <- merge(x=rat_paralogs_merged, y=rat_right,
                               by.x=c("paralog_id","paralog_name"),
                               by.y=c("paralog_id","paralog_name"), all.x=FALSE, all.y=FALSE)

rat_paralogs_merged <- rat_paralogs_merged[!duplicated(rat_paralogs_merged), ]

write.table(rat_paralogs_merged, file = "Paralogues_mapping_ppb_expression_organs-rat.txt", sep = "\t", row.names = FALSE, quote = FALSE )

rat_long <- gather(rat_paralogs_merged, condition, expression, Brain_gene:Testis_paralog, factor_key=TRUE)
rat_long$Organ <- gsub("_.*","",rat_long$condition, perl=TRUE)
rat_long$Type <- gsub(".*_","",rat_long$condition, perl=TRUE)
rat_long$Species <- rep("Rat",nrow(rat_long))
rat_long <- subset(rat_long, select=-c(condition))

rat_long_gene <- rat_long[rat_long$Type == "gene",]
rat_long_paralog <- rat_long[rat_long$Type == "paralog",]

colnames(rat_long_gene)[which(names(rat_long_gene) == "expression")] <- "Gene"
colnames(rat_long_paralog)[which(names(rat_long_paralog) == "expression")] <- "Paralog"

rat_long_all <- merge(x=rat_long_gene, y=rat_long_paralog,
                        by.x=c("paralog_id","paralog_name","gene_id","gene_name","Organ","Species"),
                        by.y=c("paralog_id","paralog_name","gene_id","gene_name","Organ","Species"),
                        all.x=FALSE, all.y=FALSE)
rat_long_all <- subset(rat_long_all, select=-c(Type.x, Type.y))
rat_long_all <- rat_long_all[order(rat_long_all$Organ),]

write.table(rat_long_all, file = "Paralogues_mapping_ppb_expression_organs-rat-longformat.txt", sep = "\t", row.names = FALSE, quote = FALSE )


human_long_all <- read.table(file = "Paralogues_mapping_ppb_expression_organs-human-longformat.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_long_all <- read.table(file = "Paralogues_mapping_ppb_expression_organs-mouse-longformat.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
rat_long_all <- read.table(file = "Paralogues_mapping_ppb_expression_organs-rat-longformat.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

human_common_organs <- subset(human_long_all, Organ == "Brain"|Organ == "Heart"|Organ == "Kidney"|Organ == "Liver"|Organ == "Lung"|Organ == "Testis")
mouse_common_organs <- subset(mouse_long_all, Organ == "Brain"|Organ == "Heart"|Organ == "Kidney"|Organ == "Liver"|Organ == "Lung"|Organ == "Testis")
rat_common_organs <- subset(rat_long_all, Organ == "Brain"|Organ == "Heart"|Organ == "Kidney"|Organ == "Liver"|Organ == "Lung"|Organ == "Testis")

common_organs_all <- rbind(human_common_organs, mouse_common_organs, rat_common_organs)


ggplot(common_organs_all, aes(x=Gene.expression..ppb., y=Paralog.expression..ppb.)) + 
  geom_bin2d(bins = 60) +
  scale_fill_continuous(type = "viridis")  + 
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Gene expression (ppb)")+
  ylab("Paralog expression (ppb)")+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
           label.y.npc = 0.2)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Paralogs")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  theme(legend.position = "none")+
  facet_grid(Species~Organ)
#facet_wrap(~Organ, nrow = 3)


