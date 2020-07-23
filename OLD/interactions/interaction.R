library(veriNA3d)
library(minipack)


#source("id_dssr_function.R")
data_dir2 <- "/orozco/projects/PDB.datamining/DSSR"


# Extract dssr from the PDB file and see nucleotide interactions
pdbID <- "1s72"
extraction <- id_dssr_extract(pdbID)
id_dssr <- extraction$id_dssr


# Extract dssr from the PDB file and see nucleotide interactions
#rna <- dssr(pdbID, exefile = "/orozco/homes/pluto/ematamoros/Programs/bin/miniconda2/bin/x3dna-dssr")

file <- paste(data_dir2, "/", toupper(pdbID), ".json", sep="")
rna <- tryCatch({
  if(file.exists(file)) {
    fromJSON(file)
    #print("reading")
  } else {
    dssr(pdbID)
  }
}, error=function(e) {
  e
})

rna_dssr <- rna$models$parameters$pairs[[1]][, c("index", "nt1", "nt2")]

# Create empty list
inter <- list()
# Loop over ids
for (i in 1:length(id_dssr)) {
  # Get id of given index "i"
  id <- id_dssr[i]
  
  # Create empty object to save interactions of given nucleotide
  interacting <- c()
  
  # Get interactions of the given nucleotide using first column of dssr data
  inds <- which(rna_dssr$nt1 == id)
  if (length(inds) > 0) {
    interacting <- append(interacting, rna_dssr[inds, "nt2"])
  }
  
  # Get interactions of the given nucleotide using second column of dssr data
  inds <- which(rna_dssr$nt2 == id)
  if (length(inds) > 0) {
    interacting <- append(interacting, rna_dssr[inds, "nt1"])
  }
  
  # Save result in list
  inter[[i]] <- interacting
}

#Create a data.frame with: seq, id_dssr, inter
table <- as.data.frame(cbind(seq = as.character(extraction$seq), id_dssr = as.character(extraction$id_dssr), inter))    #Try printout of table$inter[845] - Will print three interactions

