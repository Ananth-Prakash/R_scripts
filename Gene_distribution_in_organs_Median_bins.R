# To get the median bin values of organs as Summary after binning tissues separately
library(RColorBrewer)
library(heatmap3)
library(pheatmap)
library(Hmisc)
setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/')

dataset1 <- read.table("Binned_expression_PXD000547_Brain-CorpusCallosum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table("Binned_expression_PXD000548_Brain-AnteriorTemporalLobe.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset3 <- read.table("Binned_expression_PXD001325_Breast.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table("Binned_expression_PXD001608_Colon.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table("Binned_expression_PXD002029_Colon.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table("Binned_expression_PXD004143_Brain-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table("Binned_expression_PXD004332_Brain-PinealGland.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table("Binned_expression_PXD005819_Brain-AnteriorPituitaryGland.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table("Binned_expression_PXD006233_Brain-MiddleTemporalLobe.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset10 <- read.table("Binned_expression_PXD006675_Heart-Aorta.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset11 <- read.table("Binned_expression_PXD006675_Heart-AorticValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset12 <- read.table("Binned_expression_PXD006675_Heart-AtrialSeptum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset13 <- read.table("Binned_expression_PXD006675_Heart-InferiorVenaCava.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset14 <- read.table("Binned_expression_PXD006675_Heart-LeftAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset15 <- read.table("Binned_expression_PXD006675_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset16 <- read.table("Binned_expression_PXD006675_Heart-MitralValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset17 <- read.table("Binned_expression_PXD006675_Heart-PulmonaryArtery.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset18 <- read.table("Binned_expression_PXD006675_Heart-PulmonaryValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset19 <- read.table("Binned_expression_PXD006675_Heart-PulmonaryVein.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset20 <- read.table("Binned_expression_PXD006675_Heart-RightAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset21 <- read.table("Binned_expression_PXD006675_Heart-RightVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset22 <- read.table("Binned_expression_PXD006675_Heart-TricuspidValve.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset23 <- read.table("Binned_expression_PXD006675_Heart-VentricularSeptum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset24 <- read.table("Binned_expression_PXD008934_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset25 <- read.table("Binned_expression_PXD010154_AdiposeTissue.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset26 <- read.table("Binned_expression_PXD010154_AdrenalGland.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset27 <- read.table("Binned_expression_PXD010154_BoneMarrow.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset28 <- read.table("Binned_expression_PXD010154_Brain-PitutaryHypophysis.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset29 <- read.table("Binned_expression_PXD010154_Brain.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset30 <- read.table("Binned_expression_PXD010154_Colon.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset31 <- read.table("Binned_expression_PXD010154_Duodenum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset32 <- read.table("Binned_expression_PXD010154_Esophagus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset33 <- read.table("Binned_expression_PXD010154_FallopianTubeOviduct.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset34 <- read.table("Binned_expression_PXD010154_GallBladder.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset35 <- read.table("Binned_expression_PXD010154_Heart.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset36 <- read.table("Binned_expression_PXD010154_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset37 <- read.table("Binned_expression_PXD010154_Liver.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset38 <- read.table("Binned_expression_PXD010154_Lung.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset39 <- read.table("Binned_expression_PXD010154_LymphNode.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset40 <- read.table("Binned_expression_PXD010154_Ovary.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset41 <- read.table("Binned_expression_PXD010154_Pancreas.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset42 <- read.table("Binned_expression_PXD010154_Placenta.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset43 <- read.table("Binned_expression_PXD010154_Prostate.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset44 <- read.table("Binned_expression_PXD010154_Rectum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset45 <- read.table("Binned_expression_PXD010154_SalivaryGland.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset46 <- read.table("Binned_expression_PXD010154_SmallIntestine.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset47 <- read.table("Binned_expression_PXD010154_SmoothMuscle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset48 <- read.table("Binned_expression_PXD010154_Spleen.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset49 <- read.table("Binned_expression_PXD010154_Stomach.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset50 <- read.table("Binned_expression_PXD010154_Testis.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset51 <- read.table("Binned_expression_PXD010154_Thyroid.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset52 <- read.table("Binned_expression_PXD010154_Tonsil.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset53 <- read.table("Binned_expression_PXD010154_UrinaryBladder.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset54 <- read.table("Binned_expression_PXD010154_UterineEndometrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset55 <- read.table("Binned_expression_PXD010154_VermiformAppendix.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset56 <- read.table("Binned_expression_PXD010271_Brain-SubstantiaNigra.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset57 <- read.table("Binned_expression_PXD010271_Liver.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset58 <- read.table("Binned_expression_PXD010271_Ovary.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset59 <- read.table("Binned_expression_PXD010271_Pancreas.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset60 <- read.table("Binned_expression_PXD012131_Brain-Amygdala.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset61 <- read.table("Binned_expression_PXD012131_Brain-CaudateNucleus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset62 <- read.table("Binned_expression_PXD012131_Brain-Cerebellum.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset63 <- read.table("Binned_expression_PXD012131_Brain-EntorhinalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset64 <- read.table("Binned_expression_PXD012131_Brain-FrontalGyrus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset65 <- read.table("Binned_expression_PXD012131_Brain-InferiorParietalLobule.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset66 <- read.table("Binned_expression_PXD012131_Brain-NeoCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset67 <- read.table("Binned_expression_PXD012131_Brain-SuperiorTemporalGyrus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset68 <- read.table("Binned_expression_PXD012131_Brain-Thalamus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset69 <- read.table("Binned_expression_PXD012131_Brain-VisualCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset70 <- read.table("Binned_expression_PXD012755_Brain-CerebellarHemisphericCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset71 <- read.table("Binned_expression_PXD012755_Brain-OccipitalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset72 <- read.table("Binned_expression_PXD015079_Brain-PrefrontalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset73 <- read.table("Binned_expression_PXD015079_VermiformAppendix.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset74 <- read.table("Binned_expression_PXD020187_UmblicalArtery.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset75 <- read.table("Binned_expression_syn21443008_Brain_UniPenn-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset76 <- read.table("Binned_expression_syn21444980_Brain_Aging-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset77 <- read.table("Binned_expression_syn3606087_Brain_BLSA-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset78 <- read.table("Binned_expression_syn4624471_Brain_BLSA-Precuneus.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset79 <- read.table("Binned_expression_syn6038797_Brain_MountSinai-FrontalPole.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset80 <- read.table("Binned_expression_syn6038852_Brain_ACT-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset81 <- read.table("Binned_expression_syn7204174_Brain_Banner-DLPFC.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset82 <- read.table("Binned_expression_syn7431984_Brain_MayoClinic-TemporalCortex.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset83 <- read.table("Binned_expression_PXD008722_Heart-LeftAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset84 <- read.table("Binned_expression_PXD008722_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset85 <- read.table("Binned_expression_PXD008722_Heart-RightAtrium.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
#dataset86 <- read.table("Binned_expression_PXD012431_Breast.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


# Keep the minimum number of genes that are common in all datsets
Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene"), by.y=c("Gene"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(dataset1, dataset2)
#merged_data <- Merge_data(merged_data, dataset3)
merged_data <- Merge_data(merged_data, dataset4)
merged_data <- Merge_data(merged_data, dataset5)
merged_data <- Merge_data(merged_data, dataset6)
merged_data <- Merge_data(merged_data, dataset7)
merged_data <- Merge_data(merged_data, dataset8)
merged_data <- Merge_data(merged_data, dataset9)
merged_data <- Merge_data(merged_data, dataset10)
merged_data <- Merge_data(merged_data, dataset11)
merged_data <- Merge_data(merged_data, dataset12)
merged_data <- Merge_data(merged_data, dataset13)
merged_data <- Merge_data(merged_data, dataset14)
merged_data <- Merge_data(merged_data, dataset15)
merged_data <- Merge_data(merged_data, dataset16)
merged_data <- Merge_data(merged_data, dataset17)
merged_data <- Merge_data(merged_data, dataset18)
merged_data <- Merge_data(merged_data, dataset19)
merged_data <- Merge_data(merged_data, dataset20)
merged_data <- Merge_data(merged_data, dataset21)
merged_data <- Merge_data(merged_data, dataset22)
merged_data <- Merge_data(merged_data, dataset23)
merged_data <- Merge_data(merged_data, dataset24)
merged_data <- Merge_data(merged_data, dataset25)
merged_data <- Merge_data(merged_data, dataset26)
merged_data <- Merge_data(merged_data, dataset27)
merged_data <- Merge_data(merged_data, dataset28)
merged_data <- Merge_data(merged_data, dataset29)
merged_data <- Merge_data(merged_data, dataset30)
merged_data <- Merge_data(merged_data, dataset31)
merged_data <- Merge_data(merged_data, dataset32)
merged_data <- Merge_data(merged_data, dataset33)
merged_data <- Merge_data(merged_data, dataset34)
merged_data <- Merge_data(merged_data, dataset35)
merged_data <- Merge_data(merged_data, dataset36)
merged_data <- Merge_data(merged_data, dataset37)
merged_data <- Merge_data(merged_data, dataset38)
merged_data <- Merge_data(merged_data, dataset39)
merged_data <- Merge_data(merged_data, dataset40)
merged_data <- Merge_data(merged_data, dataset41)
merged_data <- Merge_data(merged_data, dataset42)
merged_data <- Merge_data(merged_data, dataset43)
merged_data <- Merge_data(merged_data, dataset44)
merged_data <- Merge_data(merged_data, dataset45)
merged_data <- Merge_data(merged_data, dataset46)
merged_data <- Merge_data(merged_data, dataset47)
merged_data <- Merge_data(merged_data, dataset48)
merged_data <- Merge_data(merged_data, dataset49)
merged_data <- Merge_data(merged_data, dataset50)
merged_data <- Merge_data(merged_data, dataset51)
merged_data <- Merge_data(merged_data, dataset52)
merged_data <- Merge_data(merged_data, dataset53)
merged_data <- Merge_data(merged_data, dataset54)
merged_data <- Merge_data(merged_data, dataset55)
merged_data <- Merge_data(merged_data, dataset56)
merged_data <- Merge_data(merged_data, dataset57)
merged_data <- Merge_data(merged_data, dataset58)
merged_data <- Merge_data(merged_data, dataset59)
merged_data <- Merge_data(merged_data, dataset60)
merged_data <- Merge_data(merged_data, dataset61)
merged_data <- Merge_data(merged_data, dataset62)
merged_data <- Merge_data(merged_data, dataset63)
merged_data <- Merge_data(merged_data, dataset64)
merged_data <- Merge_data(merged_data, dataset65)
merged_data <- Merge_data(merged_data, dataset66)
merged_data <- Merge_data(merged_data, dataset67)
merged_data <- Merge_data(merged_data, dataset68)
merged_data <- Merge_data(merged_data, dataset69)
merged_data <- Merge_data(merged_data, dataset70)
merged_data <- Merge_data(merged_data, dataset71)
merged_data <- Merge_data(merged_data, dataset72)
merged_data <- Merge_data(merged_data, dataset73)
merged_data <- Merge_data(merged_data, dataset74)
merged_data <- Merge_data(merged_data, dataset75)
merged_data <- Merge_data(merged_data, dataset76)
merged_data <- Merge_data(merged_data, dataset77)
merged_data <- Merge_data(merged_data, dataset78)
merged_data <- Merge_data(merged_data, dataset79)
merged_data <- Merge_data(merged_data, dataset80)
merged_data <- Merge_data(merged_data, dataset81)
merged_data <- Merge_data(merged_data, dataset82)
#merged_data <- Merge_data(merged_data, dataset83)
#merged_data <- Merge_data(merged_data, dataset84)
#merged_data <- Merge_data(merged_data, dataset85)
#merged_data <- Merge_data(merged_data, dataset86)



colnames(merged_data) <- gsub("breast", "Breast", colnames(merged_data), perl = TRUE)

write.table(merged_data, file = "Binned_intensities_all_samples.txt", sep = "\t", row.names = FALSE, quote = FALSE )


binned_samples <- merged_data

colnames(binned_samples) <- gsub("PXD010154.Sample1.PitutaryHypophysis", "PXD010154.Sample1a.PitutaryHypophysis", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.Amygdala", "PXD012131.Sample1a.Amygdala", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.CaudateNucleus", "PXD012131.Sample1b.CaudateNucleus", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.Cerebellum", "PXD012131.Sample1c.Cerebellum", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.EntorhinalCortex", "PXD012131.Sample1e.EntorhinalCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.FrontalGyrus", "PXD012131.Sample1d.FrontalGyrus", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.InferiorParietalLobule", "PXD012131.Sample1f.InferiorParietalLobule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample2.NeoCortex", "PXD012131.Sample2a.NeoCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample3.NeoCortex", "PXD012131.Sample3a.NeoCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample4.NeoCortex", "PXD012131.Sample4a.NeoCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.SuperiorTemporalGyrus", "PXD012131.Sample1g.SuperiorTemporalGyrus", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.Thalamus", "PXD012131.Sample1h.Thalamus", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012131.Sample1.VisualCortex", "PXD012131.Sample1i.VisualCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012755.Sample1.OccipitalCortex", "PXD012755.Sample1a.OccipitalCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012755.Sample2.OccipitalCortex", "PXD012755.Sample1b.OccipitalCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012755.Sample3.OccipitalCortex", "PXD012755.Sample1c.OccipitalCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012755.Sample4.OccipitalCortex", "PXD012755.Sample1d.OccipitalCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012755.Sample5.OccipitalCortex", "PXD012755.Sample1e.OccipitalCortex", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD012755.Sample6.OccipitalCortex", "PXD012755.Sample1f.OccipitalCortex", colnames(binned_samples))

colnames(binned_samples) <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                                  "Brain", colnames(binned_samples))

colnames(binned_samples) <- gsub("PXD006675.Sample1.Aorta", "PXD006675.Sample1a.Aorta", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.Aorta", "PXD006675.Sample2a.Aorta", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.Aorta", "PXD006675.Sample3a.Aorta", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.AorticValve", "PXD006675.Sample1b.AorticValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.AorticValve", "PXD006675.Sample2b.AorticValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.AorticValve", "PXD006675.Sample3b.AorticValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.AtrialSeptum", "PXD006675.Sample1c.AtrialSeptum", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.AtrialSeptum", "PXD006675.Sample2c.AtrialSeptum", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.AtrialSeptum", "PXD006675.Sample3c.AtrialSeptum", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.InferiorVenaCava", "PXD006675.Sample1d.InferiorVenaCava", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.InferiorVenaCava", "PXD006675.Sample2d.InferiorVenaCava", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.InferiorVenaCava", "PXD006675.Sample3d.InferiorVenaCava", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.LeftAtrium", "PXD006675.Sample1e.LeftAtrium", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.LeftAtrium", "PXD006675.Sample2e.LeftAtrium", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.LeftAtrium", "PXD006675.Sample3e.LeftAtrium", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.LeftVentricle", "PXD006675.Sample1f.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.LeftVentricle", "PXD006675.Sample2f.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.LeftVentricle", "PXD006675.Sample3f.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.MitralValve", "PXD006675.Sample1g.MitralValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.MitralValve", "PXD006675.Sample2g.MitralValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.MitralValve", "PXD006675.Sample3g.MitralValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.PulmonaryArtery", "PXD006675.Sample1h.PulmonaryArtery", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.PulmonaryArtery", "PXD006675.Sample2h.PulmonaryArtery", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.PulmonaryArtery", "PXD006675.Sample3h.PulmonaryArtery", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.PulmonaryValve", "PXD006675.Sample1i.PulmonaryValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.PulmonaryValve", "PXD006675.Sample2i.PulmonaryValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.PulmonaryValve", "PXD006675.Sample3i.PulmonaryValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.PulmonaryVein", "PXD006675.Sample1j.PulmonaryVein", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.PulmonaryVein", "PXD006675.Sample2j.PulmonaryVein", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.PulmonaryVein", "PXD006675.Sample3j.PulmonaryVein", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.RightAtrium", "PXD006675.Sample1k.RightAtrium", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.RightAtrium", "PXD006675.Sample2k.RightAtrium", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.RightAtrium", "PXD006675.Sample3k.RightAtrium", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.RightVentricle", "PXD006675.Sample1l.RightVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.RightVentricle", "PXD006675.Sample2l.RightVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.RightVentricle", "PXD006675.Sample3l.RightVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.TricuspidValve", "PXD006675.Sample1m.TricuspidValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.TricuspidValve", "PXD006675.Sample2m.TricuspidValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.TricuspidValve", "PXD006675.Sample3m.TricuspidValve", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample1.VentricularSeptum", "PXD006675.Sample1n.VentricularSeptum", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample2.VentricularSeptum", "PXD006675.Sample2n.VentricularSeptum", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD006675.Sample3.VentricularSeptum", "PXD006675.Sample3n.VentricularSeptum", colnames(binned_samples))


colnames(binned_samples) <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                 "Heart", colnames(binned_samples))

#cor_results <- cor(binned_samples[,-c(1)], use = "complete.obs", method = "pearson")

cor_results1 <- rcorr(as.matrix(binned_samples[,-c(1)]), type = c("pearson"))
correl_matrix <- cor_results1$r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correl_summary <- flattenCorrMatrix(cor_results1$r, cor_results1$P) 
  
Pairwise_only_Brain <- with(correl_summary, correl_summary[ grepl("Brain", row) & grepl("Brain", column),])
Pairwise_only_Brain$R2 <- Pairwise_only_Brain$cor*Pairwise_only_Brain$cor
Pairwise_only_Heart <- with(correl_summary, correl_summary[ grepl("Heart", row) & grepl("Heart", column),])
Pairwise_only_Heart$R2 <- Pairwise_only_Heart$cor*Pairwise_only_Heart$cor
Significant <- correl_summary[correl_summary$p < 0.05,]

Organs <- data.frame(ID= factor(gsub(".*\\.", "", colnames(binned_samples)[-1], perl=TRUE)))

Organ_colours <- list(ID = c(Brain = "#CB9E77", Breast="#ED9268", Colon="#ED8F6F", Heart="#A49AC8", 
                              AdiposeTissue="#66C2A5", AdrenalGland="#87B695", BoneMarrow="#A9AA86", 
                              Duodenum="#D49387", Esophagus="#BB989E", FallopianTubeOviduct="#A29CB6", 
                              GallBladder="#8F9FCA", Kidney="#B895C7", Liver="#CC90C5", Lung="#E18BC3", 
                              LymphNode="#DC96B1", Ovary="#CDA898", Pancreas="#BFB97E", Placenta="#B0CB65",
                              Prostate="#ABD851", Rectum="#BFD849", SalivaryGland="#D3D840", SmallIntestine="#E8D838", 
                              SmoothMuscle="#FCD830", Spleen="#F9D442", Stomach="#F4D059", Testis="#EECB70",
                              Thyroid="#E8C686", Tonsil="#E0C297", UrinaryBladder="#C9BAA4", UterineEndometrium="#BEB6AC",
                              VermiformAppendix="#B3B3B3", UmblicalArtery="#D4BE9E"))

annotation <- data.frame(Organs = gsub(".*\\.", "", colnames(correl_matrix), perl=TRUE))
rownames(annotation) <- colnames(correl_matrix)
annotation$Datasets <- gsub("\\..*", "", rownames(annotation), perl=TRUE)

Organ_newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation$Organs))))
Organ_annoCol <- Organ_newCols(length(unique(annotation$Organs)))
names(Organ_annoCol) <- sort(unique(annotation$Organs))
Organ_annoCol <- list(Organs = Organ_annoCol)

Dataset_newCols <- colorRampPalette(grDevices::hcl.colors(length(unique(annotation$Datasets)), palette = "Spectral"))
Dataset_annoCol <- Dataset_newCols(length(unique(annotation$Datasets)))
names(Dataset_annoCol) <- unique(annotation$Datasets)
Dataset_annoCol <- list(Datasets = Dataset_annoCol)

annoCol = c(Organ_annoCol, Dataset_annoCol)

#original publication figure
#pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
#         annotation_col = annotation, 
#         #annotation_colors = Sample_colours[1],
#         #annotation_colors = annoCol,
#         #annotation = annotation, 
#         fontsize = 8)

pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
         annotation_col = annotation, 
         annotation_colors = annoCol,
         #annotation = annotation, 
         fontsize = 8)

heatmapplot_data <- pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
         annotation_col = annotation, 
         annotation_colors = annoCol,
         #annotation = annotation, 
         fontsize = 8)

class(heatmapplot_data)
names(heatmapplot_data)

par(mar=c(2, 0, 0, 40))
heatmapplot_data$tree_row %>%
  as.dendrogram() %>%
  plot(horiz = TRUE)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

set.seed(123)
mat = matrix(rnorm(100), 10)
rownames(mat) = paste0("R", 1:10)
colnames(mat) = paste0("C", 1:10)
ha = HeatmapAnnotation(Datasets = anno_simple(1:24, pch = 1:24))
#column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
column_ha = HeatmapAnnotation(Datasets = annotation$Datasets, Organs = annotation$Organs)
#Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)
Heatmap(correl_matrix, name = " ", show_row_names = FALSE, show_column_names = FALSE, top_annotation = column_ha)




merged_data_long <- gather(merged_data, Datasets, Bins, colnames(merged_data[2]):colnames(merged_data[ncol(merged_data)]), factor_key=TRUE)
merged_data_long$Tissues <- gsub(".*\\.", "", merged_data_long$Datasets)

merged_data_long$Organs <- merged_data_long$Tissues
merged_data_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                                "Brain", merged_data_long$Organs)
merged_data_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                "Heart", merged_data_long$Organs)

#merged_data_long_aggregate <- aggregate(merged_data_long[,-c(2,4)], list("Gene" = merged_data_long$Gene, "Organs" = merged_data_long$Organs), median, na.rm =TRUE)
merged_data_long_aggregate <- aggregate(merged_data_long[,-c(2,4)], by=list(merged_data_long$Gene, merged_data_long$Organs), median, na.rm =TRUE)

merged_data_long_aggregate <- merged_data_long_aggregate[,-c(3,5)]
colnames(merged_data_long_aggregate) <- c("Gene", "Tissues", "Median_bins")

All_tissues_median_bins<- spread(merged_data_long_aggregate, Tissues, Median_bins)

colnames(All_tissues_median_bins[,-c(1)]) <- gsub("^(.*)", "Median_bins_", colnames(All_tissues_median_bins[,-c(1)]), perl=TRUE)

Gene_info <- read.table("Gene_distribution_in_organs-GeneNames-Median_intensities.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
Gene_info <- Gene_info[,c(1,2)]

All_tissues_median_bins <- merge(x=Gene_info, y=All_tissues_median_bins,
                                 by.x=c("GeneName"), by.y=c("Gene"))

# This output file is further modified to add UniProt IDs. see script Suppl1_add_UniProt_IDs.R
write.table(All_tissues_median_bins, file = "Gene_distribution_in_organs-GeneNames-Median_bin_values.txt", sep = "\t", row.names = FALSE, quote = FALSE )
