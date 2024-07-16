# To get numbers for summary table for acetylation datasets
library(dplyr)

#### Kerry's Example dataset: setwd("/Users/ananth/Documents/TPP/Plasmodium/PXD000070/FDR_0.01")

setwd("/Users/ananth/Documents/TPP/Acetylation/Comet_Batch_Scripts_and_Parameters_Files/PXD003530/")

#########################################################
###### For searches with acetylation of Alanine #########
##################### DECOY #############################
#########################################################

decoy_FDR  <- read.csv( file="Alanine/FDR_0.01/FDR_output.csv")

  decoy_FDR_0.01_PSM_count <- nrow(decoy_FDR[decoy_FDR$q_value <= 0.01,])

  decoy_FDR_0.01_PSM_unique_ModPeptides_count <- nrow(unique(decoy_FDR[decoy_FDR$q_value <= 0.01,"Peptide_mod", drop=FALSE]))

  decoy_FDR_0.01_AcetylatedPeptides <- decoy_FDR[decoy_FDR$q_value <= 0.01, c("PTM","Peptide_mod")]
  decoy_FDR_0.01_AcetylatedPeptides <- decoy_FDR_0.01_AcetylatedPeptides[grep("Acetyl", decoy_FDR_0.01_AcetylatedPeptides$PTM), ]
  decoy_FDR_0.01_AcetylatedPeptides_count <- nrow(decoy_FDR_0.01_AcetylatedPeptides)

  decoy_FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count <- nrow(unique(decoy_FDR_0.01_AcetylatedPeptides[,"Peptide_mod", drop=FALSE]))

  ###########
decoy_FLR <- read.csv( file="Alanine/FDR_0.01/Site-based_FLR_pA.csv")

  ####### 1 Model FLR method #######
  ####### Using only PTM_final_prob_q_value
  decoy_FLR_0.01_sites <- decoy_FLR[decoy_FLR$PTM_final_prob_q_value <= 0.01,]
  decoy_FLR_0.01_sites_count <- nrow(decoy_FLR_0.01_sites)
  decoy_FLR_0.01_site_unique_ModPeptides_count <- nrow(unique(decoy_FLR_0.01_sites[,"Peptide_mod", drop=FALSE]))

  decoy_FLR_0.05_sites <- decoy_FLR[decoy_FLR$PTM_final_prob_q_value <= 0.05,]
  decoy_FLR_0.05_sites_count <- nrow(decoy_FLR_0.05_sites)
  decoy_FLR_0.05_site_unique_ModPeptides_count <- nrow(unique(decoy_FLR_0.05_sites[,"Peptide_mod", drop=FALSE]))
  
  decoy_FLR_0.1_sites <- decoy_FLR[decoy_FLR$PTM_final_prob_q_value <= 0.1,]
  decoy_FLR_0.1_sites_count <- nrow(decoy_FLR_0.1_sites)
  decoy_FLR_0.1_site_unique_ModPeptides_count <- nrow(unique(decoy_FLR_0.1_sites[,"Peptide_mod", drop=FALSE]))
  
  
  ##########
decoy_Binomial <- read.csv( file="Alanine/FDR_0.01/binomial_collapsed_FLR.csv")

  decoy_Binomial_0.01_sites_count <- nrow(decoy_Binomial[decoy_Binomial$Binomial_final_prob_q_value <= 0.01,])
  decoy_Binomial_0.05_sites_count <- nrow(decoy_Binomial[decoy_Binomial$Binomial_final_prob_q_value <= 0.05,])
  decoy_Binomial_0.1_sites_count <-  nrow(decoy_Binomial[decoy_Binomial$Binomial_final_prob_q_value <= 0.1,])
  
decoy_table_MODEL_FLR_method <- data.frame(Type="Lysine + Alanine", Method="MODEL FLR method")
decoy_table_MODEL_FLR_method  <- cbind(decoy_table_MODEL_FLR_method , decoy_FDR_0.01_PSM_count,decoy_FDR_0.01_PSM_unique_ModPeptides_count, decoy_FDR_0.01_AcetylatedPeptides_count, 
                    decoy_FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count, decoy_FLR_0.01_sites_count, decoy_FLR_0.01_site_unique_ModPeptides_count,
                    decoy_FLR_0.05_sites_count, decoy_FLR_0.05_site_unique_ModPeptides_count, decoy_FLR_0.1_sites_count, decoy_FLR_0.1_site_unique_ModPeptides_count,
                    decoy_Binomial_0.01_sites_count, decoy_Binomial_0.05_sites_count, decoy_Binomial_0.1_sites_count)

colnames(decoy_table_MODEL_FLR_method) <- c("Type","Method", "FDR_0.01_PSM_count", "FDR_0.01_PSM_unique_ModPeptides_count", "FDR_0.01_AcetylatedPeptides_count", 
                             "FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count", "FLR_0.01_sites_count", "FLR_0.01_site_unique_ModPeptides_count",
                             "FLR_0.05_sites_count", "FLR_0.05_site_unique_ModPeptides_count", "FLR_0.1_sites_count", "FLR_0.1_site_unique_ModPeptides_count",
                             "Binomial_0.01_sites_count", "Binomial_0.05_sites_count", "Binomial_0.1_sites_count")


  ####### 2 Decoy FLR method ####### 
  ####### Using only q_value
  decoy_FLR_0.01_sites <- decoy_FLR[decoy_FLR$pA_q_value <= 0.01,]
  decoy_FLR_0.01_sites_count <- nrow(decoy_FLR_0.01_sites)
  decoy_FLR_0.01_site_unique_ModPeptides_count <- nrow(unique(decoy_FLR_0.01_sites[,"Peptide_mod", drop=FALSE]))

  decoy_FLR_0.05_sites <- decoy_FLR[decoy_FLR$pA_q_value <= 0.05,]
  decoy_FLR_0.05_sites_count <- nrow(decoy_FLR_0.05_sites)
  decoy_FLR_0.05_site_unique_ModPeptides_count <- nrow(unique(decoy_FLR_0.05_sites[,"Peptide_mod", drop=FALSE]))

  decoy_FLR_0.1_sites <- decoy_FLR[decoy_FLR$pA_q_value <= 0.1,]
  decoy_FLR_0.1_sites_count <- nrow(decoy_FLR_0.1_sites)
  decoy_FLR_0.1_site_unique_ModPeptides_count <- nrow(unique(decoy_FLR_0.1_sites[,"Peptide_mod", drop=FALSE]))
  
  ##########
  decoy_Binomial_0.01_sites_count <- nrow(decoy_Binomial[decoy_Binomial$A_q_value <= 0.01,])
  decoy_Binomial_0.05_sites_count <- nrow(decoy_Binomial[decoy_Binomial$A_q_value <= 0.05,])
  decoy_Binomial_0.1_sites_count <-  nrow(decoy_Binomial[decoy_Binomial$A_q_value <= 0.1,])
  
  decoy_table_DECOY_FLR_method <- data.frame(Type="Lysine + Alanine", Method="DECOY FLR method")
  decoy_table_DECOY_FLR_method  <- cbind(decoy_table_DECOY_FLR_method , decoy_FDR_0.01_PSM_count,decoy_FDR_0.01_PSM_unique_ModPeptides_count, decoy_FDR_0.01_AcetylatedPeptides_count, 
                                         decoy_FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count, decoy_FLR_0.01_sites_count, decoy_FLR_0.01_site_unique_ModPeptides_count,
                                         decoy_FLR_0.05_sites_count, decoy_FLR_0.05_site_unique_ModPeptides_count, decoy_FLR_0.1_sites_count, decoy_FLR_0.1_site_unique_ModPeptides_count,
                                         decoy_Binomial_0.01_sites_count, decoy_Binomial_0.05_sites_count, decoy_Binomial_0.1_sites_count)
  
  colnames(decoy_table_DECOY_FLR_method) <- c("Type","Method", "FDR_0.01_PSM_count", "FDR_0.01_PSM_unique_ModPeptides_count", "FDR_0.01_AcetylatedPeptides_count", 
                                              "FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count", "FLR_0.01_sites_count", "FLR_0.01_site_unique_ModPeptides_count",
                                              "FLR_0.05_sites_count", "FLR_0.05_site_unique_ModPeptides_count", "FLR_0.1_sites_count", "FLR_0.1_site_unique_ModPeptides_count",
                                              "Binomial_0.01_sites_count", "Binomial_0.05_sites_count", "Binomial_0.1_sites_count")
  

  #########################################################
  ###### For searches WITHOUT acetylation of Alanine ######
  ##################### MODEL #############################
  #########################################################

model_FDR  <- read.csv( file="WithoutAlanine/FDR_0.01/FDR_output.csv")
  
  model_FDR_0.01_PSM_count <- nrow(model_FDR[model_FDR$q_value <= 0.01,])
  
  model_FDR_0.01_PSM_unique_ModPeptides_count <- nrow(unique(model_FDR[model_FDR$q_value <= 0.01,"Peptide_mod", drop=FALSE]))
  
  model_FDR_0.01_AcetylatedPeptides <- model_FDR[model_FDR$q_value <= 0.01, c("PTM","Peptide_mod")]
  model_FDR_0.01_AcetylatedPeptides <- model_FDR_0.01_AcetylatedPeptides[grep("Acetyl", model_FDR_0.01_AcetylatedPeptides$PTM), ]
  model_FDR_0.01_AcetylatedPeptides_count <- nrow(model_FDR_0.01_AcetylatedPeptides)
  
  model_FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count <- nrow(unique(model_FDR_0.01_AcetylatedPeptides[,"Peptide_mod", drop=FALSE]))
  
  ###########
  ####### Model FLR method #######
  ####### Using only PTM_final_prob_q_value
model_FLR_sites <- read.csv( file="WithoutAlanine/FDR_0.01/Site-based_FLR_pA.csv")
  
  model_FLR_0.01_sites <- model_FLR_sites[model_FLR_sites$PTM_final_prob_q_value <= 0.01,]
  model_FLR_0.01_sites_count <- nrow(model_FLR_0.01_sites)
  model_FLR_0.01_site_unique_ModPeptides_count <- nrow(unique(model_FLR_0.01_sites[,"Peptide_mod", drop=FALSE]))
  
  model_FLR_0.05_sites <- model_FLR_sites[model_FLR_sites$PTM_final_prob_q_value <= 0.05,]
  model_FLR_0.05_sites_count <- nrow(model_FLR_0.05_sites)
  model_FLR_0.05_site_unique_ModPeptides_count <- nrow(unique(model_FLR_0.05_sites[,"Peptide_mod", drop=FALSE]))
  
  model_FLR_0.1_sites <- model_FLR_sites[model_FLR_sites$PTM_final_prob_q_value <= 0.1,]
  model_FLR_0.1_sites_count <- nrow(model_FLR_0.1_sites)
  model_FLR_0.1_site_unique_ModPeptides_count <- nrow(unique(model_FLR_0.1_sites[,"Peptide_mod", drop=FALSE]))
  
  ##########
model_Binomial <- read.csv( file="WithoutAlanine/FDR_0.01/binomial_collapsed_FLR.csv")
  
  model_Binomial_0.01_sites_count <- nrow(model_Binomial[model_Binomial$Binomial_final_prob_q_value <= 0.01,])
  model_Binomial_0.05_sites_count <- nrow(model_Binomial[model_Binomial$Binomial_final_prob_q_value <= 0.05,])
  model_Binomial_0.1_sites_count <-  nrow(model_Binomial[model_Binomial$Binomial_final_prob_q_value <= 0.1,])
  
model_combined <- data.frame(Type="Lysine", Method = "MODEL FLR method")
model_combined <- cbind(model_combined, model_FDR_0.01_PSM_count,model_FDR_0.01_PSM_unique_ModPeptides_count, model_FDR_0.01_AcetylatedPeptides_count, 
                    model_FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count, model_FLR_0.01_sites_count, model_FLR_0.01_site_unique_ModPeptides_count,
                    model_FLR_0.05_sites_count, model_FLR_0.05_site_unique_ModPeptides_count, model_FLR_0.1_sites_count, model_FLR_0.1_site_unique_ModPeptides_count,
                    model_Binomial_0.01_sites_count, model_Binomial_0.05_sites_count, model_Binomial_0.1_sites_count)

colnames(model_combined) <- c("Type", "Method", "FDR_0.01_PSM_count", "FDR_0.01_PSM_unique_ModPeptides_count", "FDR_0.01_AcetylatedPeptides_count", 
                              "FDR_0.01_AcetylatedPeptides_unique_ModPeptides_count", "FLR_0.01_sites_count", "FLR_0.01_site_unique_ModPeptides_count",
                              "FLR_0.05_sites_count", "FLR_0.05_site_unique_ModPeptides_count", "FLR_0.1_sites_count", "FLR_0.1_site_unique_ModPeptides_count",
                              "Binomial_0.01_sites_count", "Binomial_0.05_sites_count", "Binomial_0.1_sites_count")

Table <- rbind(model_combined,decoy_table_MODEL_FLR_method, decoy_table_DECOY_FLR_method)
