#!/usr/bin/Rscript

library(minipack)
library(veriNA3d)
library(bio3d)

#Load all .txt files and read on them

#path <- "/Users/eric1999/Desktop/Gitlab/trinucleotides/Mining/rmsdpdbs/04.torsionals/"
path <- "./"
all_files <- list.files(path, pattern="repres")

#Load the first ntinfo
## No em funciona Eric, disenyo estrategia alternativa.
#load(paste(path, "ntinfo_1.RData", sep=""))

path_mod <- "/orozco/projects/PDB.datamining/rmsd_data3/"
ntinfo <- list()
counter <- 1
for(i in 1:length(all_files)){
  init <- strsplit(all_files[i], split="representatives")
  init <- init[[1]][1]
  print(init)

  #Check all data in the file
  file <- readLines(paste(path, all_files[i], sep=""))
  
  
  #Check path at local computer
  path_mode <- paste(path_mod, init, "/", sep="")
  
  for(m in 1:length(file)){
    pdb_file <- read.pdb(paste(path_mode, file[m],".pdb", sep=""))
    add_nt <- measureNuc(pdb_file)
    add_nt$ntID <- file[m]
    add_nt <- add_nt[2,]
    #ntinfo <- rbind(ntinfo,add_nt)
    ntinfo[[counter]] <- add_nt
    counter <- counter + 1
  }
}
# Save output as table
ntinfo <- do.call(rbind, ntinfo)
#save(ntinfo, file = paste(path, "ntinfo_complete.RData", sep=""))
write.table(ntinfo, file=paste(path, "ntinfo.txt", sep=""), col.names=T, row.names=F, quote=F)
