##############################################################################
## Function to fit a 1D vector of continuous values into a series of gaussian
## as decided by BIC criteria
fit_data <- function(xx, num=1:4) {
    suppressPackageStartupMessages(require(mclust))
#    require(mclust)

    if (length(xx) < 2 || sd(xx) == 0) {
        return(NULL)
    }

    bic1 <- mclustBIC(xx, G=num, modelNames='V')
    gaussians <- which(max(bic1[1:length(bic1)], na.rm=TRUE) ==
                        bic1[1:length(bic1)])
    mclust1 <- mclustModel(xx, BICvalues=bic1)
    #gaussians<-mclust1$G
    pars<-mclust1$parameters
    prop<-as.matrix(pars$pro) #proportions of each gaussian
    mean<-pars$mean
    sd<-sqrt(pars$variance$sigmasq)

    seq <- seq(0, 5, length=100000)
    max_dens <- max(unlist(lapply(seq, FUN=function(x)
                                sum(dnorm(x,mean,sd)*prop))))


    referencedata <- data.frame(matrix(ncol = 4, nrow = 0))

    lapply(1:gaussians, FUN=function(index) {
            referencedata <- get(x="referencedata", envir=parent.frame(n=2))
            referencedata <- rbind(referencedata,
                                    c(prop=prop[index], mean=mean[index],
                                        sd=sd[index], max_dens=max_dens),
                                    stringsAsFactors=FALSE)
            assign(x="referencedata", value=referencedata,
                    envir=parent.frame(n=2))
    })

    colnames(referencedata) <- c("prop", "mean", "sd", "max_dens")

    return(referencedata)

}

## Function to plot the gaussians
plot_gaussians <-
function(xx=NULL, referencedata=NULL, xlim=c(1.8, 2.5), new=TRUE, ...) {

    if (is.null(referencedata)) {
        referencedata <- fit_data(xx)
    }

    if (new & !is.null(xx)) {
        hist(xx, xlim=xlim, ...)
    }
    lapply(1:nrow(referencedata),
            function(i) {
                mean <- referencedata[i, "mean"]
                sd <- referencedata[i, "sd"]
                #prop <- referencedata[i, "prop"]
                x <- seq(xlim[1], xlim[2], length=1000)
                curve(dnorm(x=x, mean=mean,
                        sd=sd), type="l", lwd=2, add=TRUE)
            })
}

## Function to plot the score obtained from the gaussians
plot_score <-
function(xx=NULL, referencedata=NULL, xlim=c(1.8, 2.5), new=TRUE, ...) {

    if (is.null(referencedata)) {
        referencedata <- fit_data(xx)
    }

    mean <- referencedata$mean
    sd <- referencedata$sd
    prop <- referencedata$prop
    max_dens <- referencedata$max_dens

    if (new & !is.null(xx)) {
        hist(xx, xlim=xlim, ...)
    }
    x <- seq(xlim[1], xlim[2], length=1000)
    scores <- unlist(lapply(x,
                        FUN=function(x) {
                                sum(dnorm(x, mean, sd) * prop/max_dens)
                        }))

    par(new = TRUE)
    plot(x, scores, col="red", lwd=4, type="l", xlab="", ylab="", axes=FALSE)
    axis(side = 4)
    mtext(text="scores", side=4, line=3, las=3, cex=3)
}

