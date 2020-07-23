#!/opt/R-3.6.0/bin/Rscript

cat("Load dependencies")
library(veriNA3d)
library(bio3d)
library(minipack)
library(parallel)

savepdb <- function(pdb_file2, eleno, filename) {
  pdb <- trim.pdb(pdb_file2, eleno=eleno)
  alts <- sort(unique(pdb$atom$alt))
  alt <- alts[!alts == c(".")][1]
  eleno <- pdb$atom$eleno[pdb$atom$alt %in% c(".", alt)]
  pdb <- trim.pdb(pdb, eleno=eleno)
  pdb$atom$eleno <- seq_len(nrow(pdb$atom))
  query3 <- paste(pdb$atom[as.character(eleno), "resno"],
                  pdb$atom[as.character(eleno), "insert"],
                  pdb$atom[as.character(eleno), "chain"],
                  sep="|")
  Unique <- unique(query3)
  resno2 <- c()
  for (i in seq_along(Unique)) {
    pdb$atom$resno[query3 == Unique[i]] <- i
  }
  pdb$atom$chain <- "A"
  pdb$atom$insert <- ""
  pdb$atom$charge <- ""
  pdb$atom$entid <- ""
  tryCatch({
    write.pdb(pdb, file=filename, segid=pdb$atom$segid)
  }, error=function(e) {
    write.pdb(pdb, file=filename, segid=pdb$atom$entid)
  })
}
cat("Dependencies loaded")

#source("/usr/people/ematamoros/trinucleotides/Neo4j_database/Extract_ntinfo_function.R")
#source("/usr/people/ematamoros/trinucleotides/R_codes/dssr_id/id_dssr_function.R")
#load("/orozco/homes/pluto/ematamoros/trinucleotides/R_codes/Trinucleotides/full_list.RData")

#All RNA PDB extract
alllist <- queryEntryList(justIDs = FALSE)
nuclist <-  alllist$pdbID[alllist$type %in% c("nuc", "prot-nuc")]
mainDir <- "/orozco/projects/PDB.datamining/rmsd_data3/"
write.table(nuclist, paste(mainDir, "nuclist.txt", sep=""), col.names=F, row.names=F, quote=F)

pdbIDs <- nuclist
pdbIDm <- nuclist
################################################
## Code to create 64 trinucleotides and their directories
bases <- c("A", "G", "C", "U")
res_unique2 <- c()
for (i in bases) {
  for (j in bases) {
    for(k in bases) {
      tri <- paste(i, j, k, sep="-")
      res_unique2 <- append(res_unique2, tri)
      newdir <- paste(mainDir, tri, sep="/")
      if (!dir.exists(newdir)) {
        cat("Creating directory: ", newdir, "\n")
        dir.create(newdir)
      }
    }
  }
}

###############################################

#Loop over all PDB IDs
#for(l in 1:length(pdbIDs)){
cores <- 32
cat("Using ", cores, " cores", "\n")
#pdbIDs2 <- pdbIDs[201:250]
mclapply(1:length(pdbIDs), mc.preschedule=FALSE, mc.cores=cores, function(l) {
  pdbID <- pdbIDs[l]
  #pdb <- pdbIDm[l]
  ###################################################################
  #Create a new column called residno
  pdb_file <-cifAsPDB(pdbID)
  if (!any(pdb_file$atom$resid %in% bases)) {
    return(NA)
  }
  print(pdbID)
  pdb_file$atom$insert[is.na(pdb_file$atom$insert)] <- "?"

  ##############################################################################
  
  #Extract ntinfo, id_dssr, localenv
  ntinfo <- extract_ntinfo(pdbID)
  id_dssr <- getID(ntinfo = ntinfo)
  
  #remove any type of break
  ind_break <- which(ntinfo$Break)
  if(length(ind_break) > 0){
    for(i in 1:length(ind_break)){
      #Case where the ind_break is found at the beginning
      if(length(ntinfo[ind_break[i] - 1, "ntID"]) == 0){
        ntinfo <- ntinfo[-(ind_break[i]):-(ind_break[i]+1),]
        id_dssr <- id_dssr[-(ind_break[i]):-(ind_break[i]+1),]
        
        #Check where the ind_break is found at the the end
      }else if(is.na(ntinfo[ind_break[i] + 1,"ntID"])){
        ntinfo <- ntinfo[-(ind_break[i]-1):-(ind_break[i]),]
        id_dssr <- id_dssr[-(ind_break[i]-1):-(ind_break[i]),]
      }else{
        ntinfo <- ntinfo[-(ind_break[i]-1):-(ind_break[i]+1),]
        id_dssr <- id_dssr[-(ind_break[i]-1):-(ind_break[i]+1),]
      }
    }
  }

  
  localenv <- ntinfo$localenv
  
  #Create dataframe with seq, id_dssr, localenv
  data <- as.data.frame(cbind(id_dssr, localenv), stringsAsFactors=F)
  

  # Add column with model number
  data$id_dssr <- as.character(data$id_dssr)
  tmpmodels <- strsplit(data$id_dssr, split=":", fixed=T)
  data$model <- unlist(lapply(tmpmodels, function(x) return(x[1])))  

  ###############################################################################

  models <- unique(data$model)
  for (model in models) {
    #Select pdb model
    pdb_file2 <- selectModel(pdb=pdb_file, model=model)

    # Because of bug in veriNA3d, reassing model to data.frame
    pdb_file2$atom$model <- model
    pdb_file2 <- veriNA3d:::.perfect_input_format(pdb_file2)

    # Create dssr ids within the pdb object
    pdb_file2$atom$id <- getID(ntinfo=pdb_file2$atom)$id_dssr

    # Get trinucleotides of the given model
    data2 <- data[data$model == model, ]
    
    #ind <- which(as.character(pdb_file2$atom$id) %in% data$id_dssr)
    
    #pdb_extraction <- pdb_file2$atom[ind,]
    #eleno <- pdb_file2$atom$eleno[ind]
    
    #a <- atom.select(pdb = pdb_file2, eleno = eleno)

    trinucs_inds <- which(data2$localenv %in% res_unique2)
    for (q in trinucs_inds) {
      t <- data2[q, "localenv"]
      rows <- unlist(lapply(q, function(x) (x-1):(x+1)))
      tri_id <- as.character(data[rows,]$id_dssr)

      #Look trinucleotides in the atom$residno previously created column
      eleno <- pdb_file2$atom$eleno[which(pdb_file2$atom$id %in% tri_id)]

      filename <- paste(mainDir, as.character(t),"/", pdbID, "_", tri_id[2], ".pdb", sep="")
      if (!file.exists(filename)) {
        ## Trim the input pdb to prepare a smaller one ------------------------
        
        savepdb(pdb_file2, eleno, filename)
      }
    }
  }
})

