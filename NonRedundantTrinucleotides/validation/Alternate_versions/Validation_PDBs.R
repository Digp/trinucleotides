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

path <- ("/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/validation/")
texts <- read_excel("/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/validation/pdb_list.xlsx")
pdb_list <- as.data.frame(texts)
pdb_list <- pdb_list[[1]]

############################Extract details for first PDB
#Try for the first PDB
pdb_pack <- minipack::validation(pdb_list[1])
pdb_pack$PDB <- pdb_list[1]

#PDB Analysis for first PDB
pdb_analysis <- c()
pdb_analysis$PDB <- pdb_list[1]
pdb_analysis$Validation_info <- "YES"
technique <- queryTechnique(pdb_list[1])
if ("Solution NMR" %in% technique){
  technique <- NA
}
resol <- queryResol(pdb_list[1])
if ("" %in% resol){
  resol <- NA
}
pdb_analysis$Technique <- technique
pdb_analysis$Resol <- resol

#Loop over all PDB ID list
for(i in 554:length(pdb_list)){
  tryCatch({
    nuc <- minipack::validation(pdb_list[i])
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
    nuc_d$PDB <- pdb_list[i]
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

  }, error=function(e){
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


texts <- read_excel("/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/validation/pdb_list.xlsx")
save(pdb_analysis, file = paste(path, "pdb_anlaysis.RData", sep=""))
save(pdb_pack, file = paste(path, "pdb_pack.RData", sep=""))

load(paste(path, "pdb_analysis_boss.RData", sep=""))
load(paste(path, "pdb_pack_boss.RData", sep=""))


pdb_validation <- rbind(pdb_pack, pdb_pack_boss)

pdb_analy <- as.data.frame(pdb_analysis)
pdb_analy_boss <- as.data.frame(pdb_analysis_boss)
pdb_analysis_total <- rbind(pdb_analy, pdb_analy_boss)

save(pdb_validation, file = paste(path, "PDB_VAL.RData", sep=""))
save(data_3, file = paste(path, "PDB_ANALYSIS.RData", sep=""))

load(paste(path, "Final_Dat.RData", sep=""))

#data <- unique(pdb_analysis_total)
#data <- rbind(data, pdb_analysis)
#data_2 <- rbind(data, pdb_analyis)

#ind <- which(grepl("NA",data_3$PDB))
#data_3 <- data_2[-ind,]
             