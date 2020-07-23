library(veriNA3d)
library(bio3d)

#Folder paths
path_folder <- "/orozco/projects/PDB.datamining/rmsd_data/"
path_subfol <- list.dirs(path_folder, full.names = TRUE, recursive = TRUE)
path_subfol <- path_subfol[-1]

#Get data working objects to RAM
data_to_RAM <- function(file, objectname) {
  assign(x=objectname, value=read.pdb(file, rm.alt=F), envir=.GlobalEnv)
}

#Loop over all folders and subfolders
for (i in 1:length(path_subfol)) {
  path_folder_x <- list.files(path_subfol[i], full.names = TRUE, recursive = TRUE)
  a <- list.files(path_subfol[i])
  len_a <- length(a)
  
  #Loop over all files in data
  invisible(mapply(path_folder_x, a, FUN=data_to_RAM))
  
  #Extract the seq_len of all the 1:len()
  row_1 <- seq_len(len_a)
  row_2 <- seq_len(len_a)
  
  #STEP 1: Get all combinations of elements of 'a', 'b', and 'c' taken 3 at a time
  temp = t(combn(c(row_1, row_2), 2))
  pdb1 <- get(a[temp[i,1]])
  pdb2 <- get(a[temp[i,2]])
  
  #Calculate RMSD
  RMSD(cif1=pdb1, cif2=pdb2)
  
}











b <- list.files("/orozco/projects/PDB.datamining/rmsd_data/A-A-A", full.names = TRUE, recursive = TRUE)
a <- length(b)
row_1 <- seq_len(a)
row_2 <- seq_len(a)


#STEP 1: Get all combinations of elements of 'a', 'b', and 'c' taken 3 at a time
temp = t(combn(c(row_1, row_2), 2))
pdb1 <- get(a[temp[i,1]])
pdb2 <- get(a[temp[i,2]])
RMSD(cif1=pdb1, cif2=pdb2)




read_files <- read.pdb("/orozco/projects/PDB.datamining/rmsd_data/A-A-A/1A4T_2:A.A13.pdb")

for (subdir in list.dirs(path_folder,recursive=FALSE)) {
  for (i in list.files(subdir)){
    print(i)
  }
}


