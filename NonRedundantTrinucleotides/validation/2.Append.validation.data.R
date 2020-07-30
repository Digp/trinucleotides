library(veriNA3d)
library(bio3d)
library(minipack)
library(dplyr)

path <-  "/orozco/projects/PDB.datamining/trinucleotides/NonRedundantTrinucleotides/validation/"
path2 <- "/orozco/projects/PDB.datamining/validations/"

#Extract the rows related with the PDB.
#REMEMBER: You need to previously update the 0.get.data document with the new PDB-ID list and info

load(paste(path, "PDB_DataType.RData", sep=""))

rna_extraction_t <- t(rna_extraction)
ind_rna <- which(rna_extraction_t[,1] >= 1)
rna_list <-rownames(rna_extraction_t[ind_rna,])

#List all files in the current path2 directory
#REMEBER: Execture the 1.obtain.updated.validations.R in order to retreive the latest PDBs validation reports
validation_files <- list.files(path2)

#####Check all files are within
validation_ids <- unlist(strsplit(validation_files, split =".txt"))
check_true <- which(which(rna_list %in% validation_ids ) == length(rna_list))

#Check validations and append data

data_complete <- vector("list", length(rna_list))

for(i in 1:length(rna_list)){
  ind <- which(grepl(rna_list[i], validation_files))
  if(length(ind) > 0){
    
    nuc_d <- c()
  
    #Extract data of interest
    nuc <- read.table(paste(path2, validation_files[ind], sep=""))
    nuc_d$id_dssr <- nuc$id_dssr
    nuc_d$suite_outlier <- nuc$pucker_outlier
    nuc_d$pucker_outlier <- nuc$suite_outlier
    nuc_d$clashes <- nuc$clashes
    nuc_d$bond_lengths <- nuc$bond_lengths
    nuc_d$bond_angles <- nuc$bond_angles
    nuc_d$chirals <- nuc$chirals
    nuc_d$planes <- nuc$planes
    nuc_d$rsrz <- nuc$rsrz
    nuc_d$PDB <- rna_list[i]
    val <- as.data.frame(nuc_d)
    
    data_complete[[i]] <- val
  }
}

#Append the data
data <- bind_rows(data_complete)

#Save the appended data
write.table(data, file = paste(path, "PDB_VALIDATION.txt", sep=""), row.names = FALSE)
