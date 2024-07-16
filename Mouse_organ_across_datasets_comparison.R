#Correlate tissue ppb expression of same organ in different datasets for mouse and rat

library(ggplot2)
library(gridExtra)

rat_heart1 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD001839/proteinGroups_ppb_final-tissue_names_PXD001839.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
rat_heart2 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD013543/proteinGroups_ppb_final-tissue_names_PXD013543.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

rat_heart1_long <- gather(rat_heart1, Samples, PXD001839_ppb, colnames(rat_heart1[3]):colnames(rat_heart1[ncol(rat_heart1)]), factor_key=TRUE)
rat_heart1_aggregate <- aggregate(rat_heart1_long[,-c(3)], by=list(rat_heart1_long$Gene.ID, rat_heart1_long$Gene.Symbol), median, na.rm =TRUE)
rat_heart1_aggregate <- rat_heart1_aggregate[,-c(3,4)]
colnames(rat_heart1_aggregate) <- c("Gene.ID","Gene.Symbol","PXD001839_Rat_Heart_ppb")

rat_heart2_long <- gather(rat_heart2, Samples, PXD013543_ppb, colnames(rat_heart2[3]):colnames(rat_heart2[ncol(rat_heart2)]), factor_key=TRUE)
rat_heart2_aggregate <- aggregate(rat_heart2_long[,-c(3)], by=list(rat_heart2_long$Gene.ID, rat_heart2_long$Gene.Symbol), median, na.rm =TRUE)
rat_heart2_aggregate <- rat_heart2_aggregate[,-c(3,4)]
colnames(rat_heart2_aggregate) <- c("Gene.ID","Gene.Symbol","PXD013543_Rat_Heart_ppb")

merged_rat_hearts <- merge(x=rat_heart1_aggregate, y=rat_heart2_aggregate,
                           by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                           all.x=FALSE, all.y=FALSE)

ggplot(merged_rat_hearts, aes(x=log2(PXD001839_Rat_Heart_ppb), y=log2(PXD013543_Rat_Heart_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Rat heart PXD001839 vs PXD013543")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

###### Mouse
## Heart
mouse_heart1 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD019394.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_heart1 <- mouse_heart1[,c(1,2,4)]
mouse_heart2 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD012636.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_heart3 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD008736.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

mouse_heart1_long <- gather(mouse_heart1, Samples, PXD019394_ppb, colnames(mouse_heart1[3]):colnames(mouse_heart1[ncol(mouse_heart1)]), factor_key=TRUE)
mouse_heart1_aggregate <- aggregate(mouse_heart1_long[,-c(3)], by=list(mouse_heart1_long$Gene.ID, mouse_heart1_long$Gene.Symbol), median, na.rm =TRUE)
mouse_heart1_aggregate <- mouse_heart1_aggregate[,-c(3,4)]
colnames(mouse_heart1_aggregate) <- c("Gene.ID","Gene.Symbol","PXD019394_Mouse_Heart_ppb")

mouse_heart2_long <- gather(mouse_heart2, Samples, PXD012636_ppb, colnames(mouse_heart2[3]):colnames(mouse_heart2[ncol(mouse_heart2)]), factor_key=TRUE)
mouse_heart2_aggregate <- aggregate(mouse_heart2_long[,-c(3)], by=list(mouse_heart2_long$Gene.ID, mouse_heart2_long$Gene.Symbol), median, na.rm =TRUE)
mouse_heart2_aggregate <- mouse_heart2_aggregate[,-c(3,4)]
colnames(mouse_heart2_aggregate) <- c("Gene.ID","Gene.Symbol","PXD012636_Mouse_Heart_ppb")

mouse_heart3_long <- gather(mouse_heart3, Samples, PXD008736_ppb, colnames(mouse_heart3[3]):colnames(mouse_heart3[ncol(mouse_heart3)]), factor_key=TRUE)
mouse_heart3_aggregate <- aggregate(mouse_heart3_long[,-c(3)], by=list(mouse_heart3_long$Gene.ID, mouse_heart3_long$Gene.Symbol), median, na.rm =TRUE)
mouse_heart3_aggregate <- mouse_heart3_aggregate[,-c(3,4)]
colnames(mouse_heart3_aggregate) <- c("Gene.ID","Gene.Symbol","PXD008736_Mouse_Heart_ppb")

merged_mouse_hearts <- merge(x=mouse_heart1_aggregate, y=mouse_heart2_aggregate,
                           by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                           all.x=FALSE, all.y=FALSE)
merged_mouse_hearts <- merge(x=merged_mouse_hearts, y=mouse_heart3_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)

#remove 0 values (Andy's suggestion)
merged_mouse_hearts <- merged_mouse_hearts[merged_mouse_hearts$PXD019394_Mouse_Heart_ppb != 0,]
merged_mouse_hearts <- merged_mouse_hearts[merged_mouse_hearts$PXD012636_Mouse_Heart_ppb != 0,]
merged_mouse_hearts <- merged_mouse_hearts[merged_mouse_hearts$PXD008736_Mouse_Heart_ppb != 0,]

mouse_heart_comp1 <- ggplot(merged_mouse_hearts, aes(x=log2(PXD019394_Mouse_Heart_ppb), y=log2(PXD012636_Mouse_Heart_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse heart PXD019394 vs PXD012636")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_heart_comp2 <- ggplot(merged_mouse_hearts, aes(x=log2(PXD019394_Mouse_Heart_ppb), y=log2(PXD008736_Mouse_Heart_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse heart PXD019394 vs PXD008736")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_heart_comp3 <- ggplot(merged_mouse_hearts, aes(x=log2(PXD012636_Mouse_Heart_ppb), y=log2(PXD008736_Mouse_Heart_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse heart PXD012636 vs PXD008736")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

grid.arrange(mouse_heart_comp1, mouse_heart_comp2, mouse_heart_comp3, ncol = 3, nrow = 1)

### Brain
mouse_brain1 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD003155.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_brain1 <- mouse_brain1[,c(1:8)]
mouse_brain2 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD005230.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_brain3 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD019394.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_brain3 <- mouse_brain3[,c(1:3)]
mouse_brain4 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD022614.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_brain5 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD004496.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


mouse_brain1_long <- gather(mouse_brain1, Samples, PXD003155_ppb, colnames(mouse_brain1[3]):colnames(mouse_brain1[ncol(mouse_brain1)]), factor_key=TRUE)
mouse_brain1_aggregate <- aggregate(mouse_brain1_long[,-c(3)], by=list(mouse_brain1_long$Gene.ID, mouse_brain1_long$Gene.Symbol), median, na.rm =TRUE)
mouse_brain1_aggregate <- mouse_brain1_aggregate[,-c(3,4)]
colnames(mouse_brain1_aggregate) <- c("Gene.ID","Gene.Symbol","PXD003155_Mouse_Brain_ppb")

mouse_brain2_long <- gather(mouse_brain2, Samples, PXD005230_ppb, colnames(mouse_brain2[3]):colnames(mouse_brain2[ncol(mouse_brain2)]), factor_key=TRUE)
mouse_brain2_aggregate <- aggregate(mouse_brain2_long[,-c(3)], by=list(mouse_brain2_long$Gene.ID, mouse_brain2_long$Gene.Symbol), median, na.rm =TRUE)
mouse_brain2_aggregate <- mouse_brain2_aggregate[,-c(3,4)]
colnames(mouse_brain2_aggregate) <- c("Gene.ID","Gene.Symbol","PXD005230_Mouse_Brain_ppb")

mouse_brain3_long <- gather(mouse_brain3, Samples, PXD019394_ppb, colnames(mouse_brain3[3]):colnames(mouse_brain3[ncol(mouse_brain3)]), factor_key=TRUE)
mouse_brain3_aggregate <- aggregate(mouse_brain3_long[,-c(3)], by=list(mouse_brain3_long$Gene.ID, mouse_brain3_long$Gene.Symbol), median, na.rm =TRUE)
mouse_brain3_aggregate <- mouse_brain3_aggregate[,-c(3,4)]
colnames(mouse_brain3_aggregate) <- c("Gene.ID","Gene.Symbol","PXD019394_Mouse_Brain_ppb")

mouse_brain4_long <- gather(mouse_brain4, Samples, PXD022614_ppb, colnames(mouse_brain4[3]):colnames(mouse_brain4[ncol(mouse_brain4)]), factor_key=TRUE)
mouse_brain4_aggregate <- aggregate(mouse_brain4_long[,-c(3)], by=list(mouse_brain4_long$Gene.ID, mouse_brain4_long$Gene.Symbol), median, na.rm =TRUE)
mouse_brain4_aggregate <- mouse_brain4_aggregate[,-c(3,4)]
colnames(mouse_brain4_aggregate) <- c("Gene.ID","Gene.Symbol","PXD022614_Mouse_Brain_ppb")

mouse_brain5_long <- gather(mouse_brain5, Samples, PXD004496_ppb, colnames(mouse_brain5[3]):colnames(mouse_brain5[ncol(mouse_brain5)]), factor_key=TRUE)
mouse_brain5_aggregate <- aggregate(mouse_brain5_long[,-c(3)], by=list(mouse_brain5_long$Gene.ID, mouse_brain5_long$Gene.Symbol), median, na.rm =TRUE)
mouse_brain5_aggregate <- mouse_brain5_aggregate[,-c(3,4)]
colnames(mouse_brain5_aggregate) <- c("Gene.ID","Gene.Symbol","PXD004496_Mouse_Brain_ppb")

merged_mouse_brains <- merge(x=mouse_brain1_aggregate, y=mouse_brain2_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)
merged_mouse_brains <- merge(x=merged_mouse_brains, y=mouse_brain3_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)
merged_mouse_brains <- merge(x=merged_mouse_brains, y=mouse_brain4_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)
merged_mouse_brains <- merge(x=merged_mouse_brains, y=mouse_brain5_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)

#remove 0 values (Andy's suggestion)
merged_mouse_brains <- merged_mouse_brains[merged_mouse_brains$PXD003155_Mouse_Brain_ppb != 0,]
merged_mouse_brains <- merged_mouse_brains[merged_mouse_brains$PXD005230_Mouse_Brain_ppb != 0,]
merged_mouse_brains <- merged_mouse_brains[merged_mouse_brains$PXD019394_Mouse_Brain_ppb != 0,]
merged_mouse_brains <- merged_mouse_brains[merged_mouse_brains$PXD022614_Mouse_Brain_ppb != 0,]
merged_mouse_brains <- merged_mouse_brains[merged_mouse_brains$PXD004496_Mouse_Brain_ppb != 0,]


mouse_brain_comp1 <- ggplot(merged_mouse_brains, aes(x=log2(PXD003155_Mouse_Brain_ppb), y=log2(PXD005230_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD003155 vs PXD005230")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp2 <- ggplot(merged_mouse_brains, aes(x=log2(PXD003155_Mouse_Brain_ppb), y=log2(PXD019394_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD003155 vs PXD019394")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp3 <- ggplot(merged_mouse_brains, aes(x=log2(PXD003155_Mouse_Brain_ppb), y=log2(PXD022614_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD003155 vs PXD022614")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp4 <- ggplot(merged_mouse_brains, aes(x=log2(PXD003155_Mouse_Brain_ppb), y=log2(PXD004496_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD003155 vs PXD004496")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp5 <- ggplot(merged_mouse_brains, aes(x=log2(PXD005230_Mouse_Brain_ppb), y=log2(PXD019394_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD005230 vs PXD019394")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp6 <- ggplot(merged_mouse_brains, aes(x=log2(PXD005230_Mouse_Brain_ppb), y=log2(PXD022614_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD005230 vs PXD022614")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp7 <- ggplot(merged_mouse_brains, aes(x=log2(PXD005230_Mouse_Brain_ppb), y=log2(PXD004496_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD005230 vs PXD004496")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp8 <- ggplot(merged_mouse_brains, aes(x=log2(PXD019394_Mouse_Brain_ppb), y=log2(PXD022614_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD019394 vs PXD022614")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp9 <- ggplot(merged_mouse_brains, aes(x=log2(PXD019394_Mouse_Brain_ppb), y=log2(PXD004496_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD019394 vs PXD004496")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_brain_comp10 <- ggplot(merged_mouse_brains, aes(x=log2(PXD022614_Mouse_Brain_ppb), y=log2(PXD004496_Mouse_Brain_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse brain PXD022614 vs PXD004496")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

grid.arrange(mouse_brain_comp1, mouse_brain_comp2, mouse_brain_comp3,
             mouse_brain_comp4, mouse_brain_comp5, mouse_brain_comp6,
             mouse_brain_comp7, mouse_brain_comp8, mouse_brain_comp9,
             mouse_brain_comp10, ncol = 3, nrow = 4)

### Liver

mouse_liver1 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD000867.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_liver2 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD003155.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_liver2 <- mouse_liver2[,c(1,2,9:14)]
mouse_liver3 <- read.table( "/Users/ananth/Documents/MaxQuant_Bechmarking/Mouse/proteinGroups_ppb_final-tissue_names_PXD019394.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
mouse_liver3 <- mouse_liver3[,c(1,2,5)]
  
mouse_liver1_long <- gather(mouse_liver1, Samples, PXD000867_ppb, colnames(mouse_liver1[3]):colnames(mouse_liver1[ncol(mouse_liver1)]), factor_key=TRUE)
mouse_liver1_aggregate <- aggregate(mouse_liver1_long[,-c(3)], by=list(mouse_liver1_long$Gene.ID, mouse_liver1_long$Gene.Symbol), median, na.rm =TRUE)
mouse_liver1_aggregate <- mouse_liver1_aggregate[,-c(3,4)]
colnames(mouse_liver1_aggregate) <- c("Gene.ID","Gene.Symbol","PXD000867_Mouse_liver_ppb")

mouse_liver2_long <- gather(mouse_liver2, Samples, PXD003155_ppb, colnames(mouse_liver2[3]):colnames(mouse_liver2[ncol(mouse_liver2)]), factor_key=TRUE)
mouse_liver2_aggregate <- aggregate(mouse_liver2_long[,-c(3)], by=list(mouse_liver2_long$Gene.ID, mouse_liver2_long$Gene.Symbol), median, na.rm =TRUE)
mouse_liver2_aggregate <- mouse_liver2_aggregate[,-c(3,4)]
colnames(mouse_liver2_aggregate) <- c("Gene.ID","Gene.Symbol","PXD003155_Mouse_liver_ppb")

mouse_liver3_long <- gather(mouse_liver3, Samples, PXD019394_ppb, colnames(mouse_liver3[3]):colnames(mouse_liver3[ncol(mouse_liver3)]), factor_key=TRUE)
mouse_liver3_aggregate <- aggregate(mouse_liver3_long[,-c(3)], by=list(mouse_liver3_long$Gene.ID, mouse_liver3_long$Gene.Symbol), median, na.rm =TRUE)
mouse_liver3_aggregate <- mouse_liver3_aggregate[,-c(3,4)]
colnames(mouse_liver3_aggregate) <- c("Gene.ID","Gene.Symbol","PXD019394_Mouse_liver_ppb")

merged_mouse_livers <- merge(x=mouse_liver1_aggregate, y=mouse_liver2_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)
merged_mouse_livers <- merge(x=merged_mouse_livers, y=mouse_liver3_aggregate,
                             by.x=c("Gene.ID","Gene.Symbol"), by.y=c("Gene.ID","Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)

#remove 0 values (Andy's suggestion)
merged_mouse_livers <- merged_mouse_livers[merged_mouse_livers$PXD000867_Mouse_liver_ppb != 0,]
merged_mouse_livers <- merged_mouse_livers[merged_mouse_livers$PXD003155_Mouse_liver_ppb != 0,]
merged_mouse_livers <- merged_mouse_livers[merged_mouse_livers$PXD019394_Mouse_liver_ppb != 0,]


mouse_liver_comp1 <- ggplot(merged_mouse_livers, aes(x=log2(PXD000867_Mouse_liver_ppb), y=log2(PXD003155_Mouse_liver_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse liver PXD000867 vs PXD003155")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_liver_comp2 <- ggplot(merged_mouse_livers, aes(x=log2(PXD000867_Mouse_liver_ppb), y=log2(PXD019394_Mouse_liver_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse liver PXD000867 vs PXD019394")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

mouse_liver_comp3 <- ggplot(merged_mouse_livers, aes(x=log2(PXD003155_Mouse_liver_ppb), y=log2(PXD019394_Mouse_liver_ppb))) + 
  geom_point(alpha=0.5) + 
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.65,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  #geom_label_repel(aes(label = Samples), size = 2, max.overlaps = 7, fill = alpha(c("white"),0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mouse liver PXD003155 vs PXD019394")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

grid.arrange(mouse_liver_comp1, mouse_liver_comp2, mouse_liver_comp3,ncol = 3, nrow = 1)

