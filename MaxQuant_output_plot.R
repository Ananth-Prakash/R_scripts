library(ggplot2)
library(reshape)

setwd("/Users/surappaa/Documents/PXD007160/")

proteingroups_final <- read.table( "ExpressionAtlas_proteinGroups_final.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(proteingroups_final) <- gsub("Reporter.intensity.corrected.", "", colnames(proteingroups_final))

# #Anterior cingulate gyrus
# #Batch-1
# colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.1", "Batch-1. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.1", "Batch-1. Alzheimer's disease - OS98-11", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.1", "Batch-1. Alzheimer's disease - OS00-11", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.1", "Batch-1. Control - A86-46", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.1", "Batch-1. Control - A93-03", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.1", "Batch-1. Parkinson's disease - OS97-53", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.1", "Batch-1. Parkinson's disease - OS02-208", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.1", "Batch-1. Alzheimer's & Parkinson's disease - OS99-30", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.1", "Batch-1. Alzheimer's & Parkinson's disease - OS00-28", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.1", "Batch-1. GIS-2", colnames(proteingroups_final))
# 
# #Batch-2
# colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.2", "Batch-2. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.2", "Batch-2. Alzheimer's disease - OS00-32", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.2", "Batch-2. Alzheimer's disease - OS03-163", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.2", "Batch-2. Control -  OS00-06", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.2", "Batch-2. Control - OS02-35", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.2", "Batch-2. Parkinson's disease - OS03-392", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.2", "Batch-2. Parkinson's disease - OS03-395", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.2", "Batch-2. Alzheimer's & Parkinson's disease - E04-49", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.2", "Batch-2. Alzheimer's & Parkinson's disease - E04-68", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.2", "Batch-2. GIS-2", colnames(proteingroups_final))
# 
# #Batch-3
# colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.3", "Batch-3. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.3", "Batch-3. Alzheimer's disease - E04-186", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.3", "Batch-3. Alzheimer's disease - E05-04", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.3", "Batch-3. Control - OS03-299", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.3", "Batch-3. Control - E05-74", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.3", "Batch-3. Parkinson's disease - OS03-396", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.3", "Batch-3. Parkinson's disease - E04-152", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.3", "Batch-3. Alzheimer's & Parkinson's disease - E05-132", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.3", "Batch-3. Alzheimer's & Parkinson's disease - E05-179", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.3", "Batch-3. GIS-2", colnames(proteingroups_final))
# 
# #Batch-4
# colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.4", "Batch-4. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.4", "Batch-4. Alzheimer's disease- E05-56", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.4", "Batch-4. Alzheimer's disease - E05-87", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.4", "Batch-4. Control - E05-130", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.4", "Batch-4. Control - E06-41", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.4", "Batch-4. Parkinson's disease - E04-169 ", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.4", "Batch-4. Parkinson's disease - E05-81", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.4", "Batch-4. Alzheimer's & Parkinson's disease - E06-15", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.4", "Batch-4. Alzheimer's & Parkinson's disease - E06-37", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.4", "Batch-4. GIS-2", colnames(proteingroups_final))
# 
# #Batch-5
# colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.5", "Batch-5. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.5", "Batch-5. Alzheimer's disease - E06-155", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.5", "Batch-5. Alzheimer's disease - E08-53", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.5", "Batch-5. Control - E08-101", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.5", "Batch-5. Control - E08-145", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.5", "Batch-5. Parkinson's disease - E09-136", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.5", "Batch-5. Parkinson's disease - E15-132", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.5", "Batch-5. Alzheimer's & Parkinson's disease - E06-128", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.5", "Batch-5. Alzheimer's & Parkinson's disease - E13-11", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.5", "Batch-5. GIS-2", colnames(proteingroups_final))
# 
# 
# #Frontal Cortex
# #Batch-1
# colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.1", "Batch-1. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.1", "Batch-1. Alzheimer's disease - E08-53", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.1", "Batch-1. Alzheimer's disease - OS00-12", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.1", "Batch-1. Control - E05-130", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.1", "Batch-1. Control - E06-41", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.1", "Batch-1. Parkinson's disease - EOS01-70", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.1", "Batch-1. Parkinson's disease - OS97-53", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.1", "Batch-1. Alzheimer's & Parkinson's disease - E06-128", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.1", "Batch-1. Alzheimer's & Parkinson's disease - E06-37", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.1", "Batch-1. GIS-2", colnames(proteingroups_final))
# 
# #Batch-2
# colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.2", "Batch-2. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.2", "Batch-2. Alzheimer's disease - OS00-32", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.2", "Batch-2. Alzheimer's disease - OS03-163", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.2", "Batch-2. Control - A86-46", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.2", "Batch-2. Control - OS03-299", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.2", "Batch-2. Parkinson's disease - OS00-18", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.2", "Batch-2. Parkinson's disease - OS97-10", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.2", "Batch-2. Alzheimer's & Parkinson's disease - OS99-30", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.2", "Batch-2. Alzheimer's & Parkinson's disease - E05-179", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.2", "Batch-2. GIS-2", colnames(proteingroups_final))
# 
# #Batch-3
# colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.3", "Batch-3. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.3", "Batch-3. Alzheimer's disease - E04-186", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.3", "Batch-3. Alzheimer's disease - E05-04", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.3", "Batch-3. Control - A93-03", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.3", "Batch-3. Control - OS02-35", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.3", "Batch-3. Parkinson's disease - E04-152", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.3", "Batch-3. Parkinson's disease - OS99-26", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.3", "Batch-3. Alzheimer's & Parkinson's disease - E05-132", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.3", "Batch-3. Alzheimer's & Parkinson's disease - OS00-28", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.3", "Batch-3. GIS-2", colnames(proteingroups_final))
# 
# #Batch-4
# colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.4", "Batch-4. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.4", "Batch-4. Alzheimer's disease - E06-155", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.4", "Batch-4. Alzheimer's disease - E05-87", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.4", "Batch-4. Control - E08-101", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.4", "Batch-4. Control - A87-50", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.4", "Batch-4. Parkinson's disease - OS02-207", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.4", "Batch-4. Parkinson's disease - OS03-396", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.4", "Batch-4. Alzheimer's & Parkinson's disease - E04-68", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.4", "Batch-4. Alzheimer's & Parkinson's disease - E06-15", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.4", "Batch-4. GIS-2", colnames(proteingroups_final))
# 
# #Batch-5
# colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.5", "Batch-5. GIS-1", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.5", "Batch-5. Alzheimer's disease - OS00-11", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.5", "Batch-5. Alzheimer's disease - OS98-11", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.5", "Batch-5. Control - OS00-06", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.5", "Batch-5. Control - E08-145", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.5", "Batch-5. Parkinson's disease - OS02-66", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.5", "Batch-5. Parkinson's disease - OS02-208", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.5", "Batch-5. Alzheimer's & Parkinson's disease - E06-13", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.5", "Batch-5. Alzheimer's & Parkinson's disease - E04-49", colnames(proteingroups_final))
# colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.5", "Batch-5. GIS-2", colnames(proteingroups_final))




#Anterior cingulate gyrus
#Batch-1
colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. GIS-2", colnames(proteingroups_final))

#Batch-2
colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. GIS-2", colnames(proteingroups_final))

#Batch-3
colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. GIS-2", colnames(proteingroups_final))

#Batch-4
colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. GIS-2", colnames(proteingroups_final))

#Batch-5
colnames(proteingroups_final) <- gsub("1.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. GIS-2", colnames(proteingroups_final))


#Frontal Cortex
#Batch-1
colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.1", "Frontal cortex_Batch-1. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.1", "Frontal cortex_Batch-1. GIS-2", colnames(proteingroups_final))

#Batch-2
colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.2", "Frontal cortex_Batch-2. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.2", "Frontal cortex_Batch-2. GIS-2", colnames(proteingroups_final))

#Batch-3
colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.3", "Frontal cortex_Batch-3. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.3", "Frontal cortex_Batch-3. GIS-2", colnames(proteingroups_final))

#Batch-4
colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.4", "Frontal cortex_Batch-4. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.4", "Frontal cortex_Batch-4. GIS-2", colnames(proteingroups_final))

#Batch-5
colnames(proteingroups_final) <- gsub("1.frontal.cortex_batch.5", "Frontal cortex_Batch-5. GIS-1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("2.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("3.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("4.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Control 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("5.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Control 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("6.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("7.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("8.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("9.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_final))
colnames(proteingroups_final) <- gsub("10.frontal.cortex_batch.5", "Frontal cortex_Batch-5. GIS-2", colnames(proteingroups_final))

proteingroups_modify <- proteingroups_final[ -c(1,2) ]
meltData <- melt(proteingroups_modify)
meltData$Batch <- meltData$variable
meltData$Batch <- gsub("\\.\\s+.*", "", meltData$Batch, perl=TRUE)
meltData$Batch <- gsub(".*_", "", meltData$Batch, perl=TRUE)
meltData$Condition <- meltData$variable
meltData$Condition <- gsub("Batch-\\d+\\.\\s", "", meltData$Condition, perl=TRUE)
#meltData$Condition <- gsub("-.*", "", meltData$Condition, perl=TRUE)
meltData$Condition <- gsub(".*_", "", meltData$Condition, perl=TRUE)
meltData$Tissue <- meltData$variable
meltData$Tissue <- gsub("_.*", "", meltData$Tissue, perl=TRUE)

ggplot( meltData, aes(x=Condition, y=value)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("PXD007160 (TMT-10plex); MaxQuant output") +
  xlab("")+
  ylab("Reporter intensity corrected")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(Tissue ~ Batch)


