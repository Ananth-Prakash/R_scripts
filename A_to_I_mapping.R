###########################################
### A to I Editing
##### UPDATED to include codon variations!
###########################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library("Biostrings")
library("dplyr")
library("stringr")
library("ggplot2")
library("data.table")

setwd('/Users/ananth/Documents/A_To_I_Editing')

max_allowed_edit_sites_per_peptide <- 1

#Read A -> I editing sites
#all_editing_sites <- read.table(file = "allSites_sorted_good.rediportal.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill =TRUE) 
all_editing_sites <- read.table(file = "/Users/ananth/Documents/A_To_I_Editing/filteredTable.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill =TRUE) 

all_editing_sites$Single_varAAcode <- all_editing_sites$varAA
all_editing_sites$Single_varAAcode <- gsub("Ala","A",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Arg","R",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Asn","N",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Asp","D",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Cys","C",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Gln","Q",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Glu","E",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Gly","G",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("His","H",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Ile","I",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Leu","L",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Lys","K",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Met","M",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Phe","F",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Pro","P",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Pyl","O",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Ser","S",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Sec","U",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Thr","T",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Trp","W",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Tyr","Y",all_editing_sites$Single_varAAcode,perl=TRUE)
all_editing_sites$Single_varAAcode <- gsub("Val","V",all_editing_sites$Single_varAAcode,perl=TRUE)

filtered_editing_sites <- all_editing_sites[all_editing_sites$Consequence != "synonymous",]
filtered_editing_sites <- filtered_editing_sites[with(filtered_editing_sites, order(UniProt,IsoformPosition)), ]

#A:check distribution of number of samples
median_nSamples <- median(filtered_editing_sites$nSamples)

ggplot(data = filtered_editing_sites, aes(x = nSamples)) +
  geom_histogram(color = "white", fill = "lightblue")+
  geom_vline(aes(xintercept = median_nSamples), color = "red", linewidth = 2)+
  scale_x_log10()+
  theme_bw()+
  ggtitle("Rediportal: A to I sites in number of samples")

#B:check the distribution of number of edit sites per protein in dataset filtered based on synonymous
synfiltered_editsite_distribution <- filtered_editing_sites %>% count(UniProt, sort = TRUE)
#C:as well as fiktered after removig evidence based on sample number cutoff (B+C)
sample_filtered_editing_sites <- filtered_editing_sites[filtered_editing_sites$nSamples >= 10,]
samplefiltered_editsite_distribution <- sample_filtered_editing_sites %>% count(UniProt, sort = TRUE)

write.table(synfiltered_editsite_distribution, file = "ForPlot1-synfiltered_editsite_distribution.txt", sep = "\t", row.names = FALSE, quote = FALSE )
write.table(samplefiltered_editsite_distribution, file = "ForPlot2-samplefiltered_editsite_distribution.txt", sep = "\t", row.names = FALSE, quote = FALSE )

##Plot1
ggplot(data = synfiltered_editsite_distribution , aes(x = n)) +
  geom_histogram(color = "white", fill = "lightgreen")+
  #geom_histogram(color = "white", fill = "lightgreen", binwidth = 1)+
  #scale_x_log10()+
  scale_y_log10()+
  xlab("Number of AtoI edit sites in each protein in REDIportal")+
  ylab("Number of proteins")+
  theme_bw()+
  ggtitle("Rediportal: number of A to I sites in each protein\n(filter: remove synonymous sites)")

##Plot2
#samplefiltered_editsite_distribution <- read.table(file = "/Users/ananth/Documents/A_To_I_Editing/1mod_per_peptide/ForPlot2-samplefiltered_editsite_distribution_additional_peptides_from_100KSet.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill =TRUE) 

ggplot(data = samplefiltered_editsite_distribution , aes(x = n)) +
  geom_histogram(color = "white", fill = "red")+
  #geom_histogram(color = "white", fill = "brown", binwidth = 1)+
  #scale_x_log10()+
  scale_y_log10()+
  xlab("Number of AtoI edit sites in each protein in REDIportal")+
  ylab("Number of proteins")+
  theme_bw()+
  ggtitle("Rediportal: number of A to I sites in each protein\n(filter: remove synonymous sites + found in samples >= 10)")

if ( !exists("FASTA_TYPE")) warning("Please specify experiment type variable: FASTA_TYPE")
#FASTA_TYPE <- "Protein"
FASTA_TYPE <- "Peptide"


#Read UniProt FASTA file
#fasta_file <- readAAStringSet("/nfs/research/juan/AtoI/Human_OneProteinPerGeneSet_May2023_UP000005640_9606.fasta")
fasta_file <- readAAStringSet("TEST.fasta")

fasta_seq <- as.data.frame(fasta_file)
fasta_seq <- tibble::rownames_to_column(fasta_seq, "name")
colnames(fasta_seq) <- c("name","seq")
fasta_seq$name <- gsub(" .*","", fasta_seq$name, perl=TRUE)

fragments <- list()

# Look for Trypsin cleavage sites ----
for(i in 1:nrow(fasta_seq)){
  
  seq_id <-fasta_seq[i,1]
  seq <- fasta_seq[i,2]
  seq_length <- nchar(fasta_seq[i,2])
  
  #get coordinates of Trypsin cleavage sites (Lysine or Arginine, but not followed by Proline)
  repl_seq <- seq
  repl_seq <- gsub("KP|RP","xx", repl_seq, perl=TRUE)
  seq_frag_coords <- data.frame(gregexpr(pattern ='K|R',repl_seq))
  
  colnames(seq_frag_coords) <- c("coord")
  #add first and last coordinates of the sequence to the list
  seq_frag_coords <- rbind(0, seq_frag_coords)
  seq_frag_coords <- rbind(seq_frag_coords, seq_length)
  
  tmp <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("peptide_id","frag", "start_cord", "end_cord","pept_len"))
  #print(seq_id)
  
  if(FASTA_TYPE == "Protein"){
    #For full-length protein
    tmp[1,"peptide_id"] <- seq_id
    tmp[1,"frag"] <- seq
    tmp[1,"start_cord"] <- 1
    tmp[1,"end_cord"] <- seq_length
    tmp[1,"pept_len"] <- seq_length
  }
  
  if(FASTA_TYPE == "Peptide"){
    #extract peptide fragments including 2 missed Trypsin cleavage sites
    # protein has at least 3 trypsin cleavage sites
    if(nrow(seq_frag_coords) > 4){
      for(j in 1:(nrow(seq_frag_coords)-3)){
        tmp[j,"peptide_id"] <- paste0(seq_id,"|Peptide_",j)
        tmp[j,"frag"] <- letter(seq, (seq_frag_coords[j,1]+1):seq_frag_coords[j+3,1])
        tmp[j,"start_cord"] <- seq_frag_coords[j,1]+1
        tmp[j,"end_cord"] <- seq_frag_coords[j+3,1]
        tmp[j,"pept_len"] <- seq_frag_coords[j+3,1]-(seq_frag_coords[j,1]+1)+1 
      }
    }
    # if protein has maximum 2 trypsin cleavage sites
    if(nrow(seq_frag_coords) <= 4){
      tmp[1,"peptide_id"] <- paste0(seq_id,"|Peptide_",1)
      tmp$frag <- seq
      tmp$start_cord <- seq_frag_coords[1,1]+1
      tmp$end_cord <- seq_frag_coords[nrow(seq_frag_coords),1]
      tmp$pept_len <- tmp$end_cord-(tmp$start_cord)+1 
    }
  }
  fragments[[i]] <- tmp
  print(paste0("Looking for Trypsin cleavage sites... ", as.character(round(i*100/nrow(fasta_seq),1)),"%"))
}

All_peptide_fragments <- do.call(rbind, fragments)
write.table(All_peptide_fragments, file = "ForPlot3-All_peptide_fragments_length_distribution.txt", sep = "\t", row.names = FALSE, quote = FALSE )

#All_peptide_fragments <- read.table(file = "/Users/ananth/Documents/A_To_I_Editing/1mod_per_peptide/ForPlot3-All_peptide_fragments_length_distribution_additional_peptides_from_100KSet.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill =TRUE) 
# Distribution of peptide fragment lengths
if(FASTA_TYPE == "Peptide"){
  median_peptlength <- median(All_peptide_fragments$pept_len)

  ##Plot3  
  ggplot(data = All_peptide_fragments , aes(x = pept_len)) +
    geom_histogram(color = "white", fill = "lightpink")+
    geom_vline(aes(xintercept = median_peptlength), color = "red", linewidth = 1)+
    geom_vline(aes(xintercept = 7), color = "blue", linewidth = 1)+
    geom_vline(aes(xintercept = 70), color = "blue", linewidth = 1)+
    scale_x_log10()+
    #scale_y_log10()+
    xlab("peptide fragment length")+
    ylab("Number of peptides")+
    theme_bw()+
    ggtitle("Peptide length distribution")
}

#All_peptide_fragments_filtered <- All_peptide_fragments[All_peptide_fragments$pept_len >=7,]
# Limit peptide fragment length to between 7 and 70 amino acids, this reduces the number of combinations on longer peptides
 All_peptide_fragments_filtered <- All_peptide_fragments[(All_peptide_fragments$pept_len >=7 & All_peptide_fragments$pept_len <= 70),]

All_peptide_fragments_filtered$UniProt <- All_peptide_fragments_filtered[,c("peptide_id")]
All_peptide_fragments_filtered$UniProt <- gsub("^tr\\||^sp\\|","",All_peptide_fragments_filtered$UniProt, perl=TRUE)
All_peptide_fragments_filtered$UniProt <- gsub("\\|.*","",All_peptide_fragments_filtered$UniProt, perl=TRUE)
All_peptide_fragments_filtered$edited  <- NA #Y or N
All_peptide_fragments_filtered$number_of_edited_sites <- NA
All_peptide_fragments_filtered$all_edit_sites <- NA
All_peptide_fragments_filtered$nSamples  <- NA
All_peptide_fragments_filtered$nTissues  <- NA
All_peptide_fragments_filtered$edited_frag <- All_peptide_fragments_filtered$frag #show edited aa in lower case?


### Match ids and coordinates between peptide fragments and AtoI sites
#### Replace fragment amino acid to AtoI edited amino acid


# For each peptide fragment
# Process peptide fragments ---- 
for(l in 1:nrow(All_peptide_fragments_filtered)){
  fragment_sites_id <- All_peptide_fragments_filtered[l,c("UniProt")]
  fragment_start_pos <- All_peptide_fragments_filtered[l,c("start_cord")]
  fragment_stop_pos  <- All_peptide_fragments_filtered[l,c("end_cord")]
  all_edited_sites <- ""
  edit_sites_subset <- filtered_editing_sites[filtered_editing_sites$UniProt == fragment_sites_id,]
  #print(fragment_sites_id)
  prev_offset <- 0
  edit_frag_list <- as.list(strsplit(All_peptide_fragments_filtered[l,c("edited_frag")],"")[[1]])
  
  if(nrow(edit_sites_subset) >= 1){
    # go through all entries of editing sites
    for(k in 1:nrow(edit_sites_subset)){
      editing_sites_id <- edit_sites_subset[k,c("UniProt")]
      editing_sites_edited_position <- edit_sites_subset[k,c("IsoformPosition")]
      editing_sites_edited_residue <- edit_sites_subset[k,c("Single_varAAcode")]
      editing_sites_nsamples <- edit_sites_subset[k,c("nSamples")]
      editing_sites_ntissues <- edit_sites_subset[k,c("nTissues")]
      
      # check if editing site falls inbetween fragment peptide
      if( (editing_sites_edited_position >= fragment_start_pos) && (editing_sites_edited_position <= fragment_stop_pos) ){
        # calculate aa replacement position on fragment peptide and replace with AtoI aa as in rediportal
        offset <- (editing_sites_edited_position - fragment_start_pos) + 1
 
        ###### NOTE:1. Some edit positions have duplicated entries in RediPortal, because of codon positions
        ###### Generate all amino acid variation for that same site
        if (offset == prev_offset){
           new_edit_residue <- paste("[",prev_edit_residue,"/",tolower(editing_sites_edited_residue),"]",sep="")
           edit_frag_list[[offset]] <- new_edit_residue
        }else{
          edit_frag_list[[offset]] <- tolower(editing_sites_edited_residue)
          prev_edit_residue <- tolower(editing_sites_edited_residue)
          prev_offset <- offset
        }
        All_peptide_fragments_filtered[l,c("edited_frag")] <- paste0(unlist(edit_frag_list),collapse="")
        all_edited_sites <- paste(all_edited_sites, editing_sites_edited_position, sep=",")
        All_peptide_fragments_filtered[l,c("all_edit_sites")] <- all_edited_sites 
        All_peptide_fragments_filtered[l,c("edited")] <- "Y"
        All_peptide_fragments_filtered[l,c("nSamples")] <- editing_sites_nsamples 
        All_peptide_fragments_filtered[l,c("nTissues")] <- editing_sites_ntissues
      }
    }
  }
  print(paste0("Processing peptide fragments... ", as.character(round(l*100/nrow(All_peptide_fragments_filtered),1)),"%"))
  
}

All_peptide_fragments_filtered$edited[is.na(All_peptide_fragments_filtered$edited)] <- "N"
All_peptide_fragments_filtered$all_edit_sites <- gsub("^,","",All_peptide_fragments_filtered$all_edit_sites)
All_peptide_fragments_filtered$number_of_edited_sites <- sapply(strsplit(All_peptide_fragments_filtered$all_edit_sites,','), uniqueN)
All_peptide_fragments_filtered$number_of_edited_sites[is.na(All_peptide_fragments_filtered$all_edit_sites)] <- 0

##################################################################################################
#### Make different versions of peptide fragments based on combinations of edited positions ######
##################################################################################################

Only_edited_fragments <- All_peptide_fragments_filtered[All_peptide_fragments_filtered$edited == "Y",]
Non_edited_fragments <- All_peptide_fragments_filtered[All_peptide_fragments_filtered$edited == "N",]

# Get only number of all possible edit sites variations in a peptide ---- 

if(FASTA_TYPE == "Peptide"){
Number_of_peptide_fragment_variations <- Only_edited_fragments
Number_of_peptide_fragment_variations$number_of_edited_sites <- sapply(strsplit(Number_of_peptide_fragment_variations$all_edit_sites,','), uniqueN)

write.table(Number_of_peptide_fragment_variations, file = "ForPlot4-Number_of_peptide_fragment_variations.txt", sep = "\t", row.names = FALSE, quote = FALSE )

#Number_of_peptide_fragment_variations <- read.table(file = "/Users/ananth/Documents/A_To_I_Editing/1mod_per_peptide/ForPlot4-Number_of_peptide_fragment_variations_additional_peptides_from_100KSet.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill =TRUE) 
##Plot4
ggplot(data = Number_of_peptide_fragment_variations, aes(x = number_of_edited_sites)) +
  geom_histogram(color = "white", fill = "brown")+
  #scale_x_log10()+
  scale_y_log10()+
  geom_vline(aes(xintercept = 20), color = "blue", linewidth = 1)+
  xlab("Number of AtoI edit sites in each peptide")+
  ylab("Number of peptides")+
  theme_bw()+
  ggtitle("Rediportal: number of A to I sites in each peptide fragment\n(filter: remove synonymous sites)")
}

# Generate edited peptide fragments versions ---- 
# different versions also include un-edited version of fragment
# ## Use expand.grid() function?
# https://stackoverflow.com/questions/18705153/generate-list-of-all-possible-combinations-of-elements-of-vector

edited_fragment_versions <- list()
generate_seq_comb <- function(seq,edited_seq){
  
  edited_seq <- data.frame(edited_seq)
  ## if edited fragment has more than one variation at an amino acid position (due to a different codon)
  ## generate variations of this peptide first
  ## Ex. peptide with 2 substitutions at the same amino acid position: ABC[d/x]EFG[h/y]IJ will result in
  ## ABCdEFGhIJ, ABCdEFGyIJ, ABCxEFGhIJ, ABCxEFGyIJ: Here each sequence variation will have 4 further variations (2^n)
  ## See NOTE1. Line 221
  seq_vars <- list()
  if(grepl("\\[", edited_seq)){
    frag_list <- as.list(strsplit(edited_seq[1,1],"\\[|\\]"))[[1]]
    frag_list_var_pos <- grep("/",frag_list)
    no_of_codon_variation_sites <- length(frag_list_var_pos)
    codon_variation_combinations <- expand.grid(rep(list(0:1), no_of_codon_variation_sites))
    for(i in 1:length(frag_list_var_pos)){
      var_residue_1 <- strsplit(frag_list[frag_list_var_pos[i]],"/")[[1]][1]
      var_residue_2 <- strsplit(frag_list[frag_list_var_pos[i]],"/")[[1]][2]
      codon_variation_combinations[codon_variation_combinations[,i]==1,i] <- var_residue_1
      codon_variation_combinations[codon_variation_combinations[,i]==0,i] <- var_residue_2
    }
    for(j in 1:nrow(codon_variation_combinations)){
      frag_list[frag_list_var_pos] <- c(t(codon_variation_combinations[j,]))
      seq_vars[[j]] <- paste0(unlist(frag_list), collapse="")
    }
    edited_seq <- do.call(rbind, seq_vars)
    #print(edited_seq)
  }
  
  for(x in 1:nrow(edited_seq)){
  
  edit_peptide_pos <- unlist(gregexpr("[[:lower:]]",edited_seq[x,]))
  no_of_edits <- length(edit_peptide_pos)
  edit_aa <- data.frame()
  output <- data.frame(edited_peptide_pos=c(NA), edited_peptide_frag=c(NA))
  
  edit_pos_combinations <- expand.grid(rep(list(0:1), no_of_edits))
  
  # 28-07-2023
  # After consulting with Juan, decided to create different versions where edits events are relaxed to 1,2,3,4,5 sites
  #
  ## FILTER for Version1 ############
  #  Assumption: Chances of a peptide having more than #### 2 RNA edits #### simultaneously is extremely unlikely invivo
  ### This reduces the computational number of peptide variations
 
  #edit_pos_combinations <- data.frame(edit_pos_combinations[rowSums(edit_pos_combinations) <= 10,])
  edit_pos_combinations <- data.frame(edit_pos_combinations[rowSums(edit_pos_combinations) <= max_allowed_edit_sites_per_peptide,])

  for(k in 1:ncol(edit_pos_combinations)){
    edit_pos_combinations[edit_pos_combinations[k] == 1, k] <- edit_peptide_pos[k]
    edit_aa[k,1] <- substr(edited_seq[x,],edit_peptide_pos[k],edit_peptide_pos[k])
  }
  
  for(n in 1:nrow(edit_pos_combinations)){
    new_seq_comb <- seq
    aa_position<-NA
    output[n,"edited_peptide_pos"] <- aa_position
    output[n,"edited_peptide_frag"] <- new_seq_comb
    for(m in 1:ncol(edit_pos_combinations)){
      if(edit_pos_combinations[n,m] != 0){
        stringr::str_sub(string = new_seq_comb, start = edit_pos_combinations[n,m], end = edit_pos_combinations[n,m]) <- edit_aa[m,1]
        aa_position <- paste(aa_position,edit_pos_combinations[n,m], sep=",")
        output[n,"edited_peptide_pos"] <- aa_position
        output[n,"edited_peptide_frag"] <- new_seq_comb
      }
    }
    output$edited_peptide_pos <- gsub("^NA,","",output$edited_peptide_pos)
   }
   edited_fragment_versions[[x]] <- output
   #print(output)
  }
   all_fragment_combinations <- do.call(rbind, edited_fragment_versions)
   all_fragment_combinations <- unique(all_fragment_combinations)
   all_fragment_combinations$var <- paste("var:",rep(seq(1:nrow(all_fragment_combinations))), sep="")
   return(all_fragment_combinations)
}

res <- list()
for(a in 1:nrow(Only_edited_fragments)){
  #seq <- "ABCDEFGHIJKL"
  #edited_seq <- "ABCD[e/x]FGHI[j/y]KL"
  seq <- Only_edited_fragments[a,c("frag")]
  edited_seq <- Only_edited_fragments[a,c("edited_frag")]

  no_of_edit_sites <- Only_edited_fragments[a,c("number_of_edited_sites")]
  
  ### FILTER: Consider only peptide that have overall upto 20 edit sites on them
  ### The more edit sites on a peptide there are more variations that will be generated 
  ### (ex. a peptide with 25 edit sites will result in 2^25 = 33,554,432 peptide variations!)
  ### See Plot4, this distribution shows only few peptides that in total have > 20 edit sites on them
  ### Therefore this limit will result in loss of only few reference peptides.
   if(no_of_edit_sites <= 20){
    frag_variations <- generate_seq_comb(seq, edited_seq)
  
    inpdata <- Only_edited_fragments[a,]
    #attach seq combinations to the orig. dataframe (by duplicating the orig. frame)
    combined <- cbind(inpdata, frag_variations)
  
    res[[a]] <- combined

    print(paste0("Generating peptide fragment variations... ", as.character(round(a*100/nrow(Only_edited_fragments),1)),"%"))
   }
}

All_edited_fragment_versions <- do.call(rbind, res)
All_edited_fragment_versions$number_of_edited_sites <- lengths(strsplit(as.character(All_edited_fragment_versions$edited_peptide_pos), ","))
All_edited_fragment_versions$number_of_edited_sites[is.na(All_edited_fragment_versions$edited_peptide_pos)] <- 0

Non_edited_fragments$"edited_peptide_pos" <- NA
Non_edited_fragments$"edited_peptide_frag" <- Non_edited_fragments$frag
Non_edited_fragments$"var" <- "var:1"

All_data <- rbind(All_edited_fragment_versions, Non_edited_fragments)

All_data <- subset(All_data, select=-c(edited_frag))
All_data$edited[is.na(All_data$edited_peptide_pos)] <- "N"

colnames(All_data)[2] <- "peptide_fragment_with_2_missed_Trp_cleavage_sites"
colnames(All_data)[9] <- "edit_position_on_protein_seq"
colnames(All_data)[12] <- "edit_position_on_peptide_fragment"
colnames(All_data)[13] <- "edited_peptide_fragment"

All_data$peptide_id <- do.call(paste, c(All_data[c("peptide_id","edit_position_on_peptide_fragment","var")], sep = "_Editpos:"))
All_data$peptide_id <- gsub("Editpos:var","var",All_data$peptide_id)

#Try to numerically order by peptide id
#ords <- data.frame(All_data$peptide_id)
#colnames(ords) <- c("char")
#ords$num <- gsub(".*_","", ords$char)
#ords$char <- gsub("_\\d+","", ords$char)
#ords$start <- All_data$start_cord
#foo <- Alldata[order(-xtfrm(ords$char), ords$num),]

Accessions <- unique(filtered_editing_sites[,c("GeneName","UniProt")])
All_data <- merge(x=Accessions, y=All_data,
                  by.x=c("UniProt"), by.y=c("UniProt"),
                  all.x=FALSE, all.y=TRUE)

All_data <- All_data[order(All_data$peptide_id),]

if(FASTA_TYPE == "Protein"){
  write.table(All_data, file = "1Mod_per_peptide_output_A_to_I_edited_FULL_Length_Proteins_TABLE.txt", sep = "\t", row.names = FALSE, quote = FALSE )} else {
    write.table(All_data, file = "1Mod_per_peptide_output_A_to_I_edited_Peptide_Fragments_TABLE.txt", sep = "\t", row.names = FALSE, quote = FALSE )
  }


output_fasta <- All_data
output_fasta$edited_peptide_fragment <- toupper(output_fasta$edited_peptide_fragment)
output_fasta$FASTA <- paste0(">",output_fasta$peptide_id,"\n",output_fasta$edited_peptide_fragment)

if(FASTA_TYPE == "Protein"){
  write.table(output_fasta$FASTA, file = "1Mod_per_peptide_output_A_to_I_edited_FULL_Length_Proteins_FASTA.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE )} else {
    write.table(output_fasta$FASTA, file = "1Mod_per_peptide_output_A_to_I_edited_Peptide_Fragments_FASTA.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE )
  }


