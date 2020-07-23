#!/opt/easybuild/software/R/3.4.4-foss-2018a-X11-20180131/bin/Rscript
#!/usr/bin/Rscript
#!/orozco/homes/pluto/ematamoros/Programs/bin/miniconda2/envs/r_env/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
    make_option(c("-p", "--path"), type="character", default=NULL,
                help="Path with PDB files to treat",
                metavar="character"),
    make_option(c("-f", "--file1"), type="character", default=NULL,
                help="PDB list",
                metavar="character"),
    make_option(c("-i", "--file2"), type="character", default=NULL,
                help="Validation",
                metavar="character"),
    make_option(c("-o", "--ofile"), type="character", default=NULL,
                help="Output file",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$path) | is.null(opt$file1) | is.null(opt$file2) | is.null(opt$ofile)){
    print_help(opt_parser)
    stop("Please, provide necessary arguments", call.=FALSE)
}

# opt <- list(path="../../../../rmsd_data3/A-A-A_out/", file1="PDB_ANALYSIS.txt", file2="PDB_VAL.txt")
## ----------------------------------------------------------------------------
## Load dependencies and helper functions
## ----------------------------------------------------------------------------
library(bio3d)
library(veriNA3d)
#library(combinat)

#Get data working objects to RAM
data_to_RAM <- function(file, objectname, pbar=NULL, i) {
    #print(file)
    assign(x=objectname, value=read.pdb(file, rm.alt=F), envir=.GlobalEnv)
#    setTxtProgressBar(pbar, i)
}

## ----------------------------------------------------------------------------
## Get data
## ----------------------------------------------------------------------------
path <- opt$path
cat(path)
pdblist <- dir(path, pattern=".pdb", full.names=T)
pdblist2 <- dir(path, pattern=".pdb", full.names=F)
pdblist2 <- gsub(".pdb", "", pdblist2)
inds <- seq(1, length(pdblist), 1)

## Read file with resolutions
resoldat <- read.table(opt$file1, header=T, stringsAsFactors=F)
inds <- grep("|", resoldat$Resol, fixed=T)
resoldat$Resol[inds] <- unlist(lapply(strsplit(resoldat[inds, "Resol"], "|", fixed=T), function(x) {return(x[1])}))

## Manage resol data and select resolutions under 2.5 AND NMR data
resoldat$Resol <- as.numeric(resoldat$Resol)
resoldat$Resol[is.na(resoldat$Resol)] <- 0
resoldat2 <- resoldat[(resoldat$Resol < 2.5), ]
resoldat2 <- resoldat2[which(resoldat2$Validation_info == "YES"),]

## Make vector of PDB IDs of the given trinucleotides
pdblist2_pdbs <- substr(pdblist2, 1, 4)
pdblist2_uni <- unique(pdblist2_pdbs)

## Filter using resolutions
pdblistA <- pdblist[pdblist2_pdbs %in% resoldat2$PDB]
pdblistA2 <- pdblist2[pdblist2_pdbs %in% resoldat2$PDB]

## Read file with validations
val <- read.table(opt$file2, header=T, stringsAsFactors=F)
## Keep data of interesting PDBs
val2 <- val[val$PDB %in% resoldat2$PDB, ]
## Generate ids with PDB IDs and ID_DSSR
val2$trinames <- paste(val2$PDB, val2$id_dssr, sep="_")

val3 <- val2[which(rowSums(val2[, c("suite_outlier", "pucker_outlier", "chirals", "clashes", "rsrz")], na.rm=T) == 0), ]
#val3 <- val2[which(val2$trinames %in% pdblistA2), ]

## Filter using validation data
inds <- which(pdblistA2 %in% val3$trinames)
pdblistB <- pdblistA[inds]
pdblistB2 <- pdblistA2[inds]

#write.table(pdblistB, paste(opt$ofile, ".list", sep=""), col.names=F, row.names=F, quote=F)
#write.table(pdblistB2, paste(opt$ofile, ".list2", sep=""), col.names=F, row.names=F, quote=F)

## Load pdbs
total <- length(pdblistB)
cat(paste("Total of ", total, " useful trinucleotides\n", sep=""))
#pbar <- txtProgressBar(min=0, max=total, style=3)
inds <- seq(1, total, 1)

#cat("Loading pdbs\n")
##invisible(mapply(pdblistB, pdblistB2, i=inds, FUN=data_to_RAM, MoreArgs=list(pbar=pbar)))
invisible(mapply(pdblistB, pdblistB2, i=inds, FUN=data_to_RAM, MoreArgs=list(pbar=NULL)))
#cat("\n")
### ----------------------------------------------------------------------------
### Calculate RMSD
### ----------------------------------------------------------------------------
combs <- t(combn(inds, 2))
#combs <- permn(inds)

total <- nrow(combs)
cat(paste("Total of ", total, " combinations\n", sep=""))
cat("Calculate RMSD\n")
#pbar <- txtProgressBar(min=0, max=total, style=3)
pbar <- txtProgressBar(min=0, max=total, style=3)

rmsdout <- c()
system.time(for (i in 1:total) {
    pdb1 <- get(pdblistB2[combs[i, 1]])
    pdb2 <- get(pdblistB2[combs[i, 2]])
    #Calculate RMSD
    rmsdout[i] <- RMSD(cif1=pdb1, cif2=pdb2)
    setTxtProgressBar(pbar, i)
})
cat("\n")

#a <- read.table(opt$ofile, header=T, stringsAsFactors=F)
#rmsdout <- a$rmsd
#combs <- t(combn(inds, 2))

out <- data.frame(pdb1=pdblistB2[combs[, 1]], pdb2=pdblistB2[combs[, 2]], rmsd=rmsdout)
write.table(out, opt$ofile, col.names=T, row.names=F, quote=F)
