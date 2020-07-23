
library(veriNA3d)
library(dplyr)
library(tidyverse)
library(readr)
library(jsonlite)

full_list <- list()
ntinfo <- list()

# Retreive all PDB IDs from the Protein Dat Bank and take only those "nuc" and "prot-nuc"

alllist <- queryEntryList(justIDs = FALSE)
nuclist <-  alllist$pdbID[alllist$type %in% c("nuc", "prot-nuc")]

#Read files from the data directory
data_dir <- "/orozco/projects/PDB.datamining/veriNA"
data_dir2 <- "/orozco/projects/PDB.datamining/DSSR"
txt_files2 <- dir(data_dir, pattern = "ntinfo.txt", full.names = TRUE)
txt_files <- dir(data_dir, pattern = "ntinfo.txt")

#Subset all PDB ID from the file names "PDB"_ntinfo.txt
pdbs <-  strsplit(txt_files, "_")
only_pdb <- list()

for (i in 1:length(pdbs)){
  only_pdb[i] <- pdbs[[i]][1] 
}

only_pdb <- toupper(as.character(only_pdb))

#-----------------------------------------------------
total <- length(nuclist)
total <- nuclist[7000:11035]
pbar <- txtProgressBar(min=0, max=total, style=3)

#Loop over all nuclist
for(i in 9500:11035){
  setTxtProgressBar(pbar, i)
  ntinfo <- tryCatch({
    if(nuclist[i] %in% only_pdb){
      a <- which(only_pdb %in% nuclist[i])
      read.table(txt_files2[as.numeric(a)], header=TRUE, stringsAsFactors=FALSE)
    } else{
      pipeNucData(nuclist[i], progressbar = FALSE)
    }
  }, error=function(e) {
    e
  })
  
  if (inherits(ntinfo, "error")) {
    next()
  }
  
  #Obtain sequence
  seq <- ntinfo$resid
  
  # Creation of resid_resno
  resid_resno <- paste(ntinfo$resid, ntinfo$resno, sep="")
  
  # Creation of chain.residresno
  chain_resid_resno <- paste(ntinfo$chain, resid_resno, sep=".")
  
  #Creation of id_dssr <- model:chain.residresno
  id_dssr <- paste(ntinfo$model, chain_resid_resno, sep=":")
  
  count <- 0
  
  for(k in ntinfo$insert){
    if(k =='?'){
      # do not do anything if that happens
      count = count + 1
    } else {
      # paste into position i of vector m
      id_dssr[count] <- paste(id_dssr[count], k, sep='^')
      count = count +1
    }
  }
  
  
  table_data <- as.data.frame(cbind(seq, id_dssr))
  
  #---------------------------------------------------------------------------
  
  motifs <- c("helix", "helix-end")
  
  #dssrdata <- dssr(nuclist[i])
  file <- paste(data_dir2, "/", toupper(nuclist[i]), ".json", sep="")
  dssrdata <- tryCatch({
    if(file.exists(file)) {
      fromJSON(file)
      #print("reading")
    } else {
      dssr(nuclist[i])
    }
  }, error=function(e) {
    e
  })
  
  if (inherits(dssrdata, "error")) {
    next()
  }
  
  pat <- list()
  data_nts_full <- list()
  count <- 1
  
  #Loop over all different models and parameters
  for(l in 1:dssrdata$num_models){
    data <- dssrdata$models$parameters$nts[[l]]
    data_nts <- data$nt_id
    summary <- data$summary
    #Extract dssr nucleotide data from each model
    data_nts_full <- append(data_nts_full, as.list(data_nts))
    
    #Loop over summary (descriptors) and save all "helix and non-helix" indexes as TRUE and others FALSE
    if (length(summary) != 0){
      for (j in 1:length(summary)){
        dat <- unlist(strsplit(summary[j], ","))
        cal <- which(dat %in% motifs)
        if (length(cal) == as.integer(0)){
          pat[count] <- list(FALSE)
        }else{
          pat[count] <- list(TRUE)
        }
        count <- count + 1
      }
     
      #Create dataframe with nucleotide data and TRUE/FALSE motif logicals
      table2 <- as.data.frame(cbind(data_nts, pat))
      
      #Create two tables with all FALSE names "motifs" and "ntinfo.localenv"
      table_data$motifs <- FALSE
      table_data$ntinfo.localenv <- FALSE
      
      #Extract those indexes of id_dssr identifiers and motif_nucleotide matches
      index <- match(table2$data_nts, table_data$id_dssr)
      
      #Length might not be the same. Check and process data
      if(length(table_data$motifs) == length(match(table2$data_nts, table_data$id_dssr))){
        table_data$motifs <- table2[match(table2$data_nts, table_data$id_dssr), ]$pat
        table_data$ntinfo.localenv <- ntinfo[match(table2$data_nts, table_data$id_dssr), ]$localenv
      }else{
        table_data$motifs[na.omit(index)] <- table2[na.omit(match(table2$data_nts, table_data$id_dssr)), ]$pat
        table_data$ntinfo.localenv[na.omit(index)] <- ntinfo[na.omit(match(table2$data_nts, table_data$id_dssr)), ]$localenv
        
        a <- as.logical(1:length(table_data$motifs))
        b <- as.logical(index)
        m <- subset(b, a)
        
        # Set as NA the nucleotides which are not found by dssr
        tryCatch({
          table_data$motifs[match(NA, m)] <- NA
          table_data$ntinfo.localenv[match(NA, m)] <- NA
        }, error=function(e) {
          e
        })
      }
      }
    }
    #Save in the list
    full_list[[i]] <- table_data
}

#Remove those which are NULL
index <- which(unlist(lapply(full_list, is.null)))
total_index <- 1:length(full_list)
total_index2 <- 1:length(full_list2)
index_not_col <- which(unlist(lapply(full_list2, ncol)) != 4)


right_index <- setdiff(total_index, index)


#Subset index from the total_list dataframe
full_list2 <- full_list[right_index]
full_list3 <- full_list2[-index_not_col]

#Bind all data
full_data <- do.call(rbind, full_list3)



