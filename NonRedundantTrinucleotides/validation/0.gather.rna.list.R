#!/opt/easybuild-0/software/R/3.4.4-goolf-1.4.10/bin/Rscript


path_lib <- "/orozco/projects/PDB.datamining/trinucleotides/NonRedundantTrinucleotides/validation/lib/"


########EXECUTED FOR THE FIRST TIME. SKIP IF YOU HAVE ALREADY
################################################################
#Obtian specific Ids for nuc or prot-nuc complexes 
id_list <- queryEntryList(justIDs = FALSE)
ind <- which(id_list$type %in% c("nuc", "prot-nuc"))
id_list <- id_list$pdbID[ind]

#Loop and extract information on the pdb
rna_extraction <- data.frame(countEntities(id_list[1]))

for(i in 2:length(id_list)){
  amb <- data.frame(countEntities(id_list[i]))
  rna_extraction <- cbind(rna_extraction, amb)
}

#Set the names to the whole created data.frame
names(rna_extraction) <- id_list

save(rna_extraction, file =paste(path_lib, "PDB_DataType.RData", sep=""))

########EXECUTED IN ORDER TO UPDATE THE PREVIOUS FILE
################################################################

#load the previously PDB_DataType to update
load(file =paste(path_lib, "PDB_DataType.RData", sep=""))

#Extract names
rna_extraction_t <- t(rna_extraction)
ready_ids <- rownames(rna_extraction_t)

#Obtian specific Ids for nuc or prot-nuc complexes 
id_list <- queryEntryList(justIDs = FALSE)
ind <- which(id_list$type %in% c("nuc", "prot-nuc"))
id_list <- id_list$pdbID[ind]

#Check the new id_list which is not found in the ready_ids
ind_matched <- which(ready_ids %in% id_list)

if((length(id_list) - length(ind_matched)) > 0){
  id_remaining <- id_list[-ind_matched]
}else{
  id_remaining <- 0
}

#Extract the new IDs details and info and append it
if(id_remaining != 0){
  for(i in 1:length(id_remaining)){
    amb <- data.frame(countEntities(id_remaining[i]))
    rna_extraction <- cbind(rna_extraction, amb)
  }
}

#Set names
if(id_remaining != 0){
  vector_names <- c(ready_ids, id_remaining)
  names(rna_extraction) <- vector_names
}

#Save the file
save(rna_extraction, file =paste(path_lib, "PDB_DataType.RData", sep=""))