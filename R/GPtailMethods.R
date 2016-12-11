##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************

print.GPtail <- function(x, ...) {
    cat("Convolution of a spline density and a Generalized Pareto\n")
    print(x$SD.x)
    cat("\n")
    cat("GPD parameters\n")
    print(x$par.y)
}


varnames <- function(object, ...) {
    UseMethod("varnames")
}

varnames.GPtail <- function(object, ...) {
    object$varnames
}

`varnames<-` <-  function(object, ..., value) {
    UseMethod("varnames<-")
}

`varnames<-.GPtail` <-  function(object, ..., value) {
    if (length(value) != 3L) {
        stop("'value' must be of length 3")
    }
    names(value) <- c("x", "y", "z")
    object$varnames <- as.character(value)
    object
}


##' Plot method for the \code{"GPtail"} class.
##'
##' @title Plot method for the \code{"GPtail"} class
##'
##' @param x An object of the \code{GPtail} class.
##'
##' @param y Not used yet.
##'
##' @param which Integer for the choice of the plot type: \code{1}
##' plots the probability densities of the thress r.vs \eqn{X}, \eqn{Y}
##' and \eqn{Z}, \code{2} plots the cumulative distribution function of
##' \eqn{Z}.
##'
##' @param ... Further arguments to be passed to \code{plot}. Not used
##' for now.
##' 
##' @return Nothing.
##'
##' @rdname plot.GPtail
##'
##' @method plot GPtail
##'
##' @S3method plot GPtail
##' 
plot.GPtail <- function(x, y = NULL, which = 1, ...) {

    fill <- TRUE
    vn <- varnames(x)
    
    ## ========================================================================
    ## Plot the densities
    ## ========================================================================
    if (which == 1) {

        if (max(sapply(vn, nchar)) < 4L) {
            labs <- paste("density of", vn, sep = " ")
        } else {
            labs <- vn
        }
        
        main <- sprintf("GPD param. for %s: sigma = %6.2f, xi = %6.2f",
                        vn[2], x$par.y["scale"], x$par.y["shape"])
        
        coll <- translude("gray", alpha = 0.9)
        
        opar <- par(mfrow = c(3, 1))
        on.exit(par(opar))
        
        par(mar = c(1.4, 4, 3, 3))
        
        plot(x$dens.x, type = "l", lwd = 2,
             main = main,
             ylab = labs[1], xlab = "",
             col = "orangered", xlim = range(x$z))

        if (fill) {
            lx <- length(x$dens.x$x)
            polygon(x = c(x$dens.x$x[1], x$dens.x$x, x$dens.x$x[lx], x$dens.x$x[lx]),
                    y = c(0, x$dens.x$y, x$dens.x$y[lx], 0),
                    col = translude("orangered", alpha = 0.5))
        }
        
        abline(h = 0, v = x$SD.x$xmax, col = coll)
        
        par(mar = c(1.4, 4, 1, 3))
        plot(x$dens.y, type = "l", col = "SteelBlue3", lwd = 2,
             ylab = labs[2], xlab = "")
        abline(h = 0, col = coll)
        
        if (fill) {
            ly <- length(x$dens.y$x)
            polygon(x = c(x$dens.y$x[1], x$dens.y$x, x$dens.y$x[ly], x$dens.y$x[ly]),
                    y = c(0, x$dens.y$y, x$dens.y$y[lx], 0),
                    col = translude("SteelBlue3", alpha = 0.5))
        }
        
        par(mar = c(3, 4, 1, 3))
      
        plot(x$dens.z, type = "l",
             ylab = labs[3], xlab = "",
             col = "SpringGreen3", lwd = 2)
        abline(h = 0, v = x$SD.x$xmax, col = coll)
        
        if (fill) {
            lz <- length(x$dens.z$x)
            polygon(x = c(x$dens.z$x[1], x$dens.z$x, x$dens.z$x[lz], x$dens.z$x[lz]),
                    y = c(0, x$dens.z$y, x$dens.z$y[lz], 0),
                    col = translude("SpringGreen3", alpha = 0.5))
        }

        ## par(opar)

    } else if (which == 2) {
  
        coll <- translude("gray", alpha = 0.9)

        main <- sprintf("GPD param. for %s: sigma = %6.2f, xi = %6.2f",
                        vn[2], x$par.y["scale"], x$par.y["shape"])
        
        plot(x$dist.z, type = "l",
             main = main,
             col = "SpringGreen3", lwd = 2, ylim = c(0, 1),
             xlab = vn[3], ylab = "cdf of Z")
        
        abline(h = c(0.0, 1.0), v = x$SD.x$xmax, col = coll)

    }  else if (which == 3) {
        
        RSLplot(data = x$ret.lev, lambda = x$lambda)

    }
    
   
    
}
