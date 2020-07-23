#!/opt/easybuild/software/R/3.4.4-foss-2018a-X11-20180131/bin/Rscript
#!/orozco/homes/pluto/ematamoros/Programs/bin/miniconda2/envs/r_env/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL,
              help="Path with PDB files to treat",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$path)){
  stop("Please, provide necessary arguments", call.=FALSE)
}

## ----------------------------------------------------------------------------
## Load dependencies
## ----------------------------------------------------------------------------
library(bio3d)
#library(veriNA3d)

#Write the functiom
checker <- function(file, outfile){
  read <- read.pdb(file, rm.alt=FALSE, verbose=FALSE)
  read$atom$eleno <- 1:nrow(read$atom)
  
  #Remove all Hydrogen atoms           
  inds <- which(grepl("[H]", read$atom$elety))
  sel <- atom.select(read, eleno=inds, inverse=TRUE)
  read <- trim.pdb(read, sel)
  #read$atom <- read$atom[-inds,]
  
  #Remove the 5' P-OH
  
  if(length(which(grepl("[P]", read$atom$elety))) >= 7){
    pinds <- which(grepl("[P]", read$atom$elety))
    inds_save <- tail(pinds, n=6)
    inds_rem <- setdiff(pinds, inds_save)
    sel <- atom.select(read, eleno=inds_rem, inverse=TRUE)
    read <- trim.pdb(read, sel)
  }
  write.pdb(read, file=outfile)
}

path <- opt$path
pdblist <- dir(path, full.names=T)
pdblist2 <- dir(path, full.names=F)

outdir <- paste(path, "_out", sep="")
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

total <- length(pdblist)
#pbar <- txtProgressBar(min=0, max=total, style=3)

for (i in 1:total) {
  file <- pdblist[i]
  outfile <- paste(outdir, "/", pdblist2[i], sep="")
  checker(file, outfile=outfile)
  #setTxtProgressBar(pbar, i)
}

cat("\n")
