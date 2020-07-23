library(minipack)
library(bio3d)
library(veriNA3d)
library(lattice)
library(rbmn)
library(BBmisc)
library(readxl)
library(xlsx)
library(readtext)

`%nin%` = Negate(`%in%`)

lpath <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/validation/"
texts <- readLines(paste(lpath, "pdblist.txt", sep=""))
pdb_list_2 <- texts


############################Extract details for first PDB
#Try for the first PDB
ind <- which(grepl(pdb_list_2[1],files_included))
nuc <- read.table(paste(path, files_included[ind], sep =""))
#pdb_pack <- minipack::validation(pdb_list[1758])
nuc_d <- c()
nuc_d$seq <- nuc$seq
nuc_d$id_dssr <- nuc$id_dssr
nuc_d$suite_outlier <- nuc$pucker_outlier
nuc_d$pucker_outlier <- nuc$suite_outlier
nuc_d$clashes <- nuc$clashes
nuc_d$bond_lengths <- nuc$bond_lengths
nuc_d$bond_angles <- nuc$bond_angles
nuc_d$chirals <- nuc$chirals
nuc_d$planes <- nuc$planes
nuc_d$rsrz <- nuc$rsrz
nuc_d$PDB <- pdb_list_2[1]
nuc_e <- as.data.frame(nuc_d)
pdb_pack <- nuc_e


#Loop over all PDB ID list
for(i in 3:length(pdb_list_2)){
  tryCatch({
    amb <- which(grepl(pdb_list_2[i],files_included))
    if(length(amb) > 0){
      ind <- which(grepl(pdb_list_2[i],files_included))
      nuc <- read.table(paste(path, files_included[ind], sep =""))
      
      nuc_d <- c()
      
      nuc_d$seq <- nuc$seq
      nuc_d$id_dssr <- nuc$id_dssr
      nuc_d$suite_outlier <- nuc$pucker_outlier
      nuc_d$pucker_outlier <- nuc$suite_outlier
      nuc_d$clashes <- nuc$clashes
      nuc_d$bond_lengths <- nuc$bond_lengths
      nuc_d$bond_angles <- nuc$bond_angles
      nuc_d$chirals <- nuc$chirals
      nuc_d$planes <- nuc$planes
      nuc_d$rsrz <- nuc$rsrz
      nuc_d$PDB <- pdb_list_2[i]
      nuc_e <- as.data.frame(nuc_d)
      
      ping <- which(names(pdb_pack) %nin% names(nuc_e))
      if(length(ping) > 0){
        miss <- names(pdb_pack[ping])
        
        cant <- c()
        
        for(i in 1:length(miss)){
          cant[i] <- NA
        }
        names(cant) <- miss
        nuc_m <- append(cant, nuc_e)
        nuc_m <- as.data.frame(nuc_m)
      }else if(length(ping) == 0){
        nuc_m <- nuc_e
      }
      
      
      
      #Bind with already existing information
      pdb_pack <- rbind(pdb_pack, nuc_m)
      
      #Put YES in the validation info
      pdb_analysis$PDB[i] <- pdb_list[i]
      pdb_analysis$Validation_info[i] <- "YES"
      technique <- queryTechnique(pdb_list[i])
      if ("Solution NMR" %in% technique){
        technique <- NA
      }
      resol <- queryResol(pdb_list[i])
      if ("" %in% resol){
        resol <- NA
      }
      pdb_analysis$Technique[i] <- technique
      pdb_analysis$Resol[i] <- resol
    }else{
      nuc <- minipack::validation(pdb_list_2[i])
  
      nuc_d <- c()
      nuc_d$seq <- nuc$seq
      nuc_d$id_dssr <- nuc$id_dssr
      nuc_d$suite_outlier <- nuc$pucker_outlier
      nuc_d$pucker_outlier <- nuc$suite_outlier
      nuc_d$clashes <- nuc$clashes
      nuc_d$bond_lengths <- nuc$bond_lengths
      nuc_d$bond_angles <- nuc$bond_angles
      nuc_d$chirals <- nuc$chirals
      nuc_d$planes <- nuc$planes
      nuc_d$rsrz <- nuc$rsrz
      nuc_d$PDB <- pdb_list_2[i]
      nuc_e <- as.data.frame(nuc_d)
      
    
      ping <- which(names(pdb_pack) %nin% names(nuc_e))
      if(length(ping) > 0){
        miss <- names(pdb_pack[ping])
        
        cant <- c()
        
        for(i in 1:length(miss)){
          cant[i] <- NA
        }
        names(cant) <- miss
        nuc_m <- append(cant, nuc_e)
        nuc_m <- as.data.frame(nuc_m)
      }else if(length(ping) == 0){
        nuc_m <- nuc_e
      }
      
      #Bind with already existing information
      pdb_pack <- rbind(pdb_pack, nuc_m)
    }
    
    #Put YES in the validation info
    pdb_analysis$PDB[i] <- pdb_list[i]
    pdb_analysis$Validation_info[i] <- "YES"
    technique <- queryTechnique(pdb_list[i])
    if ("Solution NMR" %in% technique){
      technique <- NA
    }
    resol <- queryResol(pdb_list[i])
    if ("" %in% resol){
      resol <- NA
    }
    pdb_analysis$Technique[i] <- technique
    pdb_analysis$Resol[i] <- resol

  }, error=function(e){
    print(paste("Error in:", pdb_list_2[i], sep=""))
   
     #Put YES in the validation info
    pdb_analysis$PDB[i] <- pdb_list[i]
    pdb_analysis$Validation_info[i] <- "NO"
    technique <- queryTechnique(pdb_list[i])
    if ("Solution NMR" %in% technique){
      technique <- NA
    }
    resol <- queryResol(pdb_list[i])
    if ("" %in% resol){
      resol <- NA
    }
    pdb_analysis$Technique[i] <- technique
    pdb_analysis$Resol[i] <- resol

  })

}

pdb_analysis_boss <- pdb_analysis
pdb_pack_boss <- pdb_pack

save(pdb_pack, file = paste(lpath, "Final_Dat.RData", sep=""))

load(paste(lpath, "PDB_VAL.RData", sep=""))
load(paste(lpath, "Anlaysis_all.RData", sep=""))

