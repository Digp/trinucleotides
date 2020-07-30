library(minipack)
library(veriNA3d)

path <- "/orozco/projects/PDB.datamining/validations/"

all_files <- list.files(path)

id_list <- queryEntryList(justIDs = FALSE)
ind <- which(id_list$type %in% c("nuc", "prot-nuc"))
id_list <- id_list$pdbID[ind]

#Extract data for every PDB

#Check if the PDBS are rna or not
for(i in 1:length(id_list)){
  if(length(which(grepl(id_list[i], all_files))) > 0){
    next()
  }else{
    tryCatch({
    print(i)
    pdb <- minipack::validation(id_list[i])
    write.table(pdb, file = paste(path, id_list[i],".txt", sep=""), row.names = FALSE)
    },
    error = function(e){
      print(paste(i, "not worked", sep=" "))
    })
  }
}


