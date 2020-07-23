#' Obtain nucleotide details from a data set of RNA structures
#' 
#' Pipeline to generate a data.frame with the desired info for a list of PDB. 
#' Nucleotides are labeled with a unique identifier (column ntID).
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#'
#' @return A data.frame with data about every nucleotide in the input set.
#'
#' @examples 
#'     ## This is a toy example, see vignettes for real-world usages.
#'     pdblist <- list("1bau", "2rn1")
#'     ntinfo <- extract_ntinfo(pdbID=pdblist)
#'
#' @author Eric Matamoros & Diego Gallego
#'

extract_ntinfo <- function(pdbID){
  #Read files from the data directory
  data_dir <- "/orozco/projects/PDB.datamining/veriNA/"
  txt_files2 <- dir(data_dir, pattern = "ntinfo.txt", full.names = TRUE)
  txt_files <- dir(data_dir, pattern = "ntinfo.txt")
  
  #Subset all PDB ID from the file names "PDB"_ntinfo.txt
  pdbs <-  strsplit(txt_files, "_")
  only_pdb <- list()
  
  for (i in 1:length(pdbs)){
    only_pdb[i] <- pdbs[[i]][1] 
  }
  
  only_pdb <- toupper(as.character(only_pdb))
  
  #-----------------------------------------------------------------------------------------------------
  #Download PDB from Proein Data Bank or retreive it from file depending on whether it is there or not
  if (toupper(pdbID) %in% only_pdb){
    a <- which(only_pdb %in% toupper(pdbID))
    ntinfo <- read.table(txt_files2[as.numeric(a)], header=TRUE, stringsAsFactors=FALSE)
  }else{
    ntinfo <- pipeNucData(pdbID = pdbID, range = c(0, 1000000))
  }
}

