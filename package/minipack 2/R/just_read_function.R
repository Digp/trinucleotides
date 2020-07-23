just_read <- function(pdbID = NULL, pdb = NULL, cif = NULL){
  if(!is.null(pdb)){
    file <- read.pdb(pdb)
  }else if(!is.null(cif)){
    file <- read.cif(cif)
    #Buggy, not working 
  }else if(!is.null(pdbID)){
    file <- pipeNucData(toupper(pdbID))
  }
  return(file)
}
