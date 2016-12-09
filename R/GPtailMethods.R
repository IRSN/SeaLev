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
    
    if (which == 1) {
        
        coll <- translude("gray", alpha = 0.9)
        
        opar <- par(mfrow = c(3, 1))
        par(mar = c(1.4, 4, 3, 3))
        
        plot(x$dens.x, type = "l", lwd = 2,
             main = sprintf("sigma_Y = %6.2f, xi_Y = %6.2f",
                 x$par.y["scale"], x$par.y["shape"]),
             ylab = "density of X", xlab = "",
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
             ylab = "density of Y", xlab = "")
        abline(h = 0, col = coll)
        if (fill) {
            ly <- length(x$dens.y$x)
            polygon(x = c(x$dens.y$x[1], x$dens.y$x, x$dens.y$x[ly], x$dens.y$x[ly]),
                    y = c(0, x$dens.y$y, x$dens.y$y[lx], 0),
                    col = translude("SteelBlue3", alpha = 0.5))
        }
        
        par(mar = c(3, 4, 1, 3))
        plot(x$dens.z, type = "l",
             ylab = "density of Z", xlab = "",
             col = "SpringGreen3", lwd = 2)
        abline(h = 0, v = x$SD.x$xmax, col = coll)
        if (fill) {
            lz <- length(x$dens.z$x)
            polygon(x = c(x$dens.z$x[1], x$dens.z$x, x$dens.z$x[lz], x$dens.z$x[lz]),
                    y = c(0, x$dens.z$y, x$dens.z$y[lz], 0),
                    col = translude("SpringGreen3", alpha = 0.5))
        }

        par(opar)

    } else if (which == 2) {
        coll <- translude("gray", alpha = 0.9)
        plot(x$dist.z, type = "l",
             main = sprintf("sigma_Y = %6.2f, xi_Y = %6.2f",
                 x$par.y["scale"], x$par.y["shape"]),
             col = "SpringGreen3", lwd = 2, ylim = c(0, 1),
             xlab = "z", ylab = "cdf of Z")
        abline(h = c(0.0, 1.0), v = x$SD.x$xmax, col = coll)
    }  else if (which == 3) {
        RSLplot(data = x$ret.lev, lambda = x$lambda)
    }
    
    


    
}
