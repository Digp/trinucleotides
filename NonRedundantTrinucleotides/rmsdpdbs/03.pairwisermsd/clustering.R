#!/opt/easybuild/software/R/3.4.4-foss-2018a-X11-20180131/bin/Rscript
#!/usr/bin/Rscript
#!/orozco/homes/pluto/ematamoros/Programs/bin/miniconda2/envs/r_env/bin/Rscript

## ----------------------------------------------------------------------------
## Read arguments
## ----------------------------------------------------------------------------
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input file with rmsd calculations",
                metavar="character"),
    make_option(c("-o", "--ofile"), type="character", default=NULL,
                help="Output file",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$ofile)){
    print_help(opt_parser)
    stop("Please, provide necessary arguments", call.=FALSE)
}

# opt <- list(input="A-A-A.txt", ofile="A-A-A.representatives.txt")
## ----------------------------------------------------------------------------
## Load dependencies and helper functions
## ----------------------------------------------------------------------------
aswide <- function(pwc){
    labels <- unique(append(unique(pwc$pdb1), unique(pwc$pdb2)))

    M <- matrix(0, nrow=length(labels), ncol=length(labels))
    colnames(M)<-labels
    rownames(M)<-labels

    for (i in 1:nrow(pwc)) {
        print(i)
        x <- pwc[i, ]
        M[which(labels==x[, 1]), which(labels==x[, 2])] <- x[, 3]
        M[which(labels==x[, 2]), which(labels==x[, 1])] <- x[, 3]
    }

    return(M)
}

pwcalcs <- read.table(opt$input, header=T, stringsAsFactors=F)
pw_wide <- aswide(pwcalcs)

h <- 0.85
hc_RMSD <- hclust(as.dist(pw_wide), method="average")
cllist_RMSD <- cutree(hc_RMSD,h=h)
clusters <- unique(cllist_RMSD)

representants <- c()
for (i in 1:length(clusters)) {
    cl <- clusters[i]
    representants[i] <- names(cllist_RMSD[which(cllist_RMSD == cl)][1])
}

write.table(representants, opt$ofile, col.names=F, row.names=F, quote=F)
