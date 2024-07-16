# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")

library(mygene)



data.to.map <- tmp
data.to.map$"ENSG" <- "NA"
data.to.map$"Gene.Symbol" <- "NA"
data.to.map$"unique.gene.count" <- "NA"


for(i in 1:nrow(data.to.map)){
  
  x <- data.frame(strsplit(data.to.map[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
  x_temp <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  f = file()
  sink(file=f)
  res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("human")),
                  error = function(e) {print(0)})
  sink()
  close(f)
  
  if (class(res)=="DFrame" | class(res) == "DataFrame"){
    data.to.map[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    if(data.to.map[i,"ENSG"] == ""){
      data.to.map[ i, "ENSG"] <- paste( unique(unlist(res$ensembl[!is.na(res$ensembl)])), collapse = ";")
    }
    data.to.map[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")}
  temp_symb <- data.to.map[i,"Gene.Symbol"]
  data.to.map[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
  
  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(data.to.map),1)),"%"))
}


#take a backup of data before filtering
data.to.map_before_filtering <- data.to.map

# remove protein groups that have no mapping to an ENSG gene IDs
data.to.map <- data.to.map[ data.to.map$ENSG != "" , ]
data.to.map <- data.to.map[ data.to.map$ENSG != "NA" , ]

# remove all protein groups that map to multiple ENSG gene IDs (this removes a lot of proteins) - 
# the reasoning to remove these cases is that we cannot establish for sure which gene is contributing the signal to the protein abundance; all genes contribute equally or one gene is a majority? 
data.to.map <- data.to.map[ data.to.map[, "unique.gene.count"] == 1, ]
data.to.map <- data.to.map[ grep(";", data.to.map$ENSG, invert = TRUE) , ]
