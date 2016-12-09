##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************

##' Plot method for \code{SplineDensity} objects.
##'
##' The density value is computed on a fine grid using the
##' \code{predict} method for the class. The original points used in
##' the creation of the object are shown as well.
##' 
##' @title Plot method for \code{SplineDensity} objects
##'
##' @param x A \code{SplineDensity} object.
##'
##' @param y Not used
##'
##' @param ... 
##'
##' @return Nothing.
##'

plot.SplineDensity <- function(x, y = NULL, ...) {

    plot(x$x, x$f, type = "l", col = "orangered", lwd = 2, ...)

    abline(v = x$knots0, col = translude("gray", alpha = 0.7))
    abline(h = 0)

    lines(x$x, x$fp, col = translude("SpringGreen3", alpha = 0.4),
          type = "o", pch = 16, cex = 0.8,
          lty = "dashed", lwd = 2)
    
    NULL
    
}

summary.SplineDensity <- function(object, ...) {
    
    class(objet) <- "summary.SplineDensity"
}

print.SplineDensity <- function(x, ...) {

    cat("SplineDensity object\n")
    cat("o order (degree + 1) ", x$order, "\n")
    cat("o Number of knots    ", length(x$knots0), "\n")
    
    
}


##' Predict (or evaluate) a Spline Density.
##'
##' The density or one of its derivatives or the cumulative
##' distribution function is evaluated at the points given in
##' \code{newdata}.
##' 
##' @title Predict (or evaluate) a Spline Density.
##'
##' @param object A \code{SplineDensity} objet.
##' 
##' @param newdata Numeric vector giving the values where the density
##' (or its derivative) will be evaluated.
##'
##' @param deriv Integer giving the a derivation order. It can be set
##' to \eqn{-1} to get the indefinite integral.
##'
##' @param ... 
##'
##' @return A list containing two numeric vectors \code{x} and
##' \code{y} representing the prediction abscissae and the
##' corresponding predicted values.
##'
##' @examples
##' set.seed(1234)
##' SD <- rSplineDensity(order = 4, xmax = 10)
##' xg <- seq(from = -1, to = 11, length.out = 200)
##' p <- predict(SD, newdata = xg)
##' plot(p, type = "l", main = "density")
##' pm1 <- predict(SD, newdata = xg, deriv = -1)
##' plot(pm1, type = "l", main = "cdf")
##' 
predict.SplineDensity <- function(object, newdata = NULL, deriv = 0, ...) {

    check <- FALSE
    
    if (is.null(newdata)) {
        xr <- object$xmax - object$xmin
        newdata <- seq(from = object$xmin - xr / 10,
                       to = object$xmax + xr / 10,
                       length.out = 200)
    } else {
        newdata <- sort(newdata)
    }
    
    f <- rep(0, length(newdata))
    ind <- (newdata > object$xmin) & (newdata < object$xmax)

    if (deriv >= 0) {
        
        ind <- (newdata > object$xmin) & (newdata < object$xmax)

        if (any(ind)) {
            x <- newdata[ind]
            B <-  splineDesign(knots = object$knots,
                               x = x,
                               ord = object$order,
                               derivs = rep(deriv, length(x)), 
                               outer.ok = TRUE)
            f[ind] <- B %*% object$b
            
        }
        
    } else if (deriv == -1) {
        
        ind <- (newdata < object$xmax)
        f[!ind] <- 1.0
        
        if (any(ind)) {
            x <- newdata[ind]
            k <- object$order
            k1 <- k - 1L
            nKnots <- length(object$knots0)
            knotsPlus <- object$knots
            knotsPlus <- c(knotsPlus, knotsPlus[length(knotsPlus)])
            
            BInt <- splineDesign(knots = knotsPlus,
                                 x = x,
                                 ord = k + 1L,
                                 outer.ok = TRUE)
            
            B  <- matrix(0.0, nrow = length(x), ncol = nKnots + k - 2L)
            dk <- rep(NA, nKnots + k1 - 1L)
            
            ##=======================================================================
            ## integration formula. See Carl de Boor (2001) 'A
            ## Practical Guide to Splines', Revised edition,
            ## Springer-Verlag, p. 128.
            ##
            ## XXX TODO Correct the computation so that integration begins at xmin
            ## rather than at the first knot which is slightly smaller.
            ## =======================================================================
            
            ## first, compute the integrals from the first knot to 'xmax'
            for (ell in 1L:(nKnots + k1 - 1L)) {
                dk[ell] <- (knotsPlus[ell + k] - knotsPlus[ell]) / k
                B[ , ell] <- dk[ell] * apply(BInt[ , ell:(nKnots + k1 - 1L), drop = FALSE], 1, sum)
            }
            
            if (check) {
                matplot(x, B, type = "l")
            }
            
            f[ind] <- B %*% object$b
        }
        
    } else {
        stop("only -1 is allowed as negative derivation order for now")
    }
            
    
    list(x = newdata, y = f)

}
