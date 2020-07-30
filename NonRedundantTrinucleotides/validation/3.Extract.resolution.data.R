library(veriNA3d)
library(bio3d)
library(minipack)
library(dplyr)

path <-  "/orozco/projects/PDB.datamining/trinucleotides/NonRedundantTrinucleotides/validation/"

########EXECUTED FOR THE FIRST TIME. SKIP IF YOU HAVE ALREADY
################################################################
load(paste(path, "PDB_DataType.RData", sep=""))

rna_extraction_t <- t(rna_extraction)
ind_rna <- which(rna_extraction_t[,1] >= 1)
rna_list <-rownames(rna_extraction_t[ind_rna,])

pdb_analysis <- c()
#Loop over the remaining PDBs
for(i in 1:length(rna_list)){
  
  #Check for technique and resolution
  technique <- queryTechnique(rna_list[i])
  if ("Solution NMR" %in% technique){
    technique <- NA
  }
  resol <- queryResol(rna_list[i])
  if ("" %in% resol || is.null(resol) ){
    resol <- NA
  }
  pdb_analysis$Technique[i] <- technique
  pdb_analysis$Resol[i] <- resol
  pdb_analysis$PDB[i] <- rna_list[i]
  pdb_analysis$Validation_info[i] <- "YES"
  
}

pdb_analysis <- as.data.frame(pdb_analysis)

#Save the data
write.table(pdb_analysis, file = paste(path, "PDB_ANALYSIS.txt", sep=""), row.names = FALSE)


########EXECUTED IN ORDER TO UPDATE THE PREVIOUS FILE
################################################################

#Load the previous existing file
read_val <- read.table(paste(path, "PDB_ANALYSIS.txt", sep=""))
names(read_val) <- c("Technique", "Resol", "PDB", "Validation_info")
read_val <- read_val[-1,]


#Extract the rows related with the PDB.
#REMEMBER: You need to previously update the 0.get.data document with the new PDB-ID list and info

load(paste(path, "PDB_DataType.RData", sep=""))

rna_extraction_t <- t(rna_extraction)
ind_rna <- which(rna_extraction_t[,1] >= 1)
rna_list <-rownames(rna_extraction_t[ind_rna,])

#Select those ids which are not found
ind <- which(rna_list %in% read_val$PDB)
rna_remaining <- rna_list[-ind]

pdb_analysis <- c()

#Loop over the remaining PDBs
for(i in 1:length(rna_remaining)){
  
  #Check for technique and resolution
  technique <- queryTechnique(rna_remaining[i])
  if ("Solution NMR" %in% technique){
    technique <- NA
  }
  resol <- queryResol(rna_remaining[i])
  if ("" %in% resol || is.null(resol) ){
    resol <- NA
  }
  pdb_analysis$Technique[i] <- technique
  pdb_analysis$Resol[i] <- resol
  pdb_analysis$PDB[i] <- rna_remaining[i]
  pdb_analysis$Validation_info[i] <- "YES"
  
}

pdb_analysis <- as.data.frame(pdb_analysis)

#Append to the original data

data_tec <- rbind(read_val, pdb_analysis)

#Save the data
write.table(data_tec, file = paste(path, "PDB_ANALYSIS.txt", sep=""), row.names = FALSE)


