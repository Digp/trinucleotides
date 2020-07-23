library(veriNA3d)
library(dplyr)
library(tidyverse)
library(readr)

alllist <- queryEntryList(justIDs = FALSE)
nuclist <-  alllist$pdbID[alllist$type %in% c("nuc", "prot-nuc")]


data_dir <- "/orozco/projects/PDB.datamining/veriNA/"
txt_files2 <- dir(data_dir, pattern = "ntinfo.txt", full.names = TRUE)
txt_files <- dir(data_dir, pattern = "ntinfo.txt")

#pdbID <- c("100D", "1BAU")
#pdbID <- as.list(pdbID)

pdbs <-  strsplit(txt_files, "_")
only_pdb <- list()

for (i in 1:length(pdbs)){
  only_pdb[i] <- pdbs[[i]][1] 
}

only_pdb <- toupper(as.character(only_pdb))
ntinfo2 <- list()
#which(only_pdb %in% nuclist)

for(i in 1:length(nuclist)){
  tryCatch({
    if(nuclist[i] %in% only_pdb){
      a <- which(only_pdb %in% nuclist[i])
      ntinfo2[[i]] <- read.table(txt_files2[as.numeric(a)], header=TRUE, stringsAsFactors=FALSE)$localenv
    } else{
      ntinfo2[[i]] <- checkNuc(pdb = nuclist[i])$localenv
    }
  }, error=function(e) {
    #next()
  })
}

#Remove all identifiers that do contain aberrant modifications or extremes
ntd <- unlist(ntinfo2)[grep("^[AUCG]-[AUGC]-[AUGC]$", unlist(ntinfo2), perl=T)]
#Nucleotide <- droplevels(ntd)
Nucleotide <- ntd
plat <- as.data.frame(table(Nucleotide))

plat2 <- plat %>%
  arrange(desc(Freq))

#Create ggplot. Seems not to sort data in a descending or ascending way
ggplot(plat2, aes(x = Nucleotide, y = Freq)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle=90,margin = margin(1, unit = "cm"),vjust =1))

#-----------------------------------

#Create the data frame with ascending order data
sorted <- data.frame(reorder(plat2$Nucleotide, plat2$Freq))
names(sorted) <- paste("Nucleotide")
sorted_data <- data.frame (Nucleotide = sorted, Freq = plat2$Freq)

#Crete the plot with ascending order data
ggplot(sorted_data, aes(x = Nucleotide, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90,margin = margin(1, unit = "cm"),vjust =1))

