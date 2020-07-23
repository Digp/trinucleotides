#' Extract the unique identifier for each of the nucleotides
#' 
#' Pipeline to generate the ID input for the dssr which contains
#' the format: model.nucleotide:chain^subchain. Result output as string
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param ntinfo A data.frame containing the information retreived from the
#'     PipeNucData function from the veriNa3d package containing unique
#'     parameters for each nucleotide
#'     
#'
#' @return A string with the unique identifier for each value. If input a list
#'     or a data.frame column the output results in a column with all single identifiers
#'
#' @examples 
#'     ## This is a toy example, see vignettes for real-world usages.
#'     pdblist <- list("1bau", "2rn1")
#'     ntinfo <- extract_ntinfo(pdbID=pdblist)
#'     id_dssr <- id_dssr_extract(ntinfo)
#'     ntinfo$id_dssr <- id_dssr   #Both id_dssr and ntinfo have the same length
#'
#' @author Eric Matamoros

#Extract PDBID data and 
id_dssr_extract <- function(pdbID=NULL, ntinfo=NULL){

  #-------------------------------------------------------
  
  if(is.null(ntinfo) ==  TRUE){
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
  data_table <- as.data.frame(cbind(seq, id_dssr))
  }
  return(data_table)
} 
  
