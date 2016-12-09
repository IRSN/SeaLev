##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************

##' Random Spline Density.
##'
##' From the sequence of knots and the order, a basis of B-splines is
##' built. A vector of positive random coefficients is drawn and
##' suitably normalized to ensure that the spline integrates to one on
##' the support.
##'
##' @title Random Spline Density
##'
##' @param knots Numeric vector of knots.
##'
##' @param nx Number of default evaluation points.
##'
##' @param nKnots Number of knots if \code{knots} is not given.
##'
##' @param xmin Lower end-point of the desnity.
##'
##' @param xmax Upper end-point of the density.
##'
##' @param order Order of the spline, i.e. the degree of the
##' polynomial pieces minus one. With \code{ord = 2} one get a broken
##' line spline and \code{ord = 4} leads to a cubic spline.
##'
##' @param plot Logical. If \code{TRUE} a plot is shown.
##'
##' @return An object of class \code{"SplineDensity"}.
##'
##' @author Yves Deville
##'
##' @examples
##' set.seed(1234)
##' SD <- rSplineDensity(order = 4, xmax = 10)
##' plot(SD, main = "spline density")
##' ## uniform knots sequence
##' SDu <- rSplineDensity(order = 4, knots = 0:10)
##' plot(SDu, main = "spline density with uniform knots sequence.")
rSplineDensity <- function(knots,
                           nx = 200L,
                           nKnots = 10L,
                           xmin = 0.0,
                           xmax = 1.0,
                           order = 4L,
                           plot = TRUE){

    if (missing(knots)) {
        xR <- xmax - xmin
        knots <- cumsum(runif(nKnots - 1))
        knots <- xmin + c(0, knots / knots[length(knots)]) * xR
    } else {
        nKnots <- length(knots)
        xmin <- min(knots)
        xmax <- max(knots)
        xR <- xmax - xmin
    }
    
    k <- order
    k1 <- k - 1L
    x <-  seq(from = xmin, to = xmax, length.out = nx)
    h <- x[2L] - x[1L]
    
    knotsPlus <- knots
    eps <- 1e-3 * xR
    knotsPlus <- c(rep(knots[1] - eps, k1), knots, rep(knots[nKnots] + eps, k1))
  
    B <- splineDesign(knots = knotsPlus,
                      x = x,
                      ord = k,
                      outer.ok = TRUE)
    
    Bd <- splineDesign(knots = knotsPlus,
                       x = knotsPlus,
                       ord = k,
                       derivs = rep(k1 , length(knotsPlus)),
                       outer.ok = TRUE)

    knotsPlusPlus <- c(knotsPlus, knotsPlus[length(knotsPlus)])
    BInt <- splineDesign(knots = knotsPlusPlus,
                         x = c(xmin, xmax),
                         ord = k + 1L,
                         outer.ok = TRUE)

    ints  <- matrix(0.0, nrow = 2L, ncol = nKnots + k - 2L)
    dk <- rep(NA, nKnots + k1 - 1L)
    
    ##=======================================================================
    ## integration formula. See
    ##
    ## Carl de Boor (2001) 'A Practical Guide to Splines',
    ## Revised edition, Springer-Verlag, p. 128.
    ##
    ##=======================================================================
    
    ## first, compute the integrals from the first knot to 'xmax'
    for (ell in 1L:(nKnots + k1 - 1L)) {
        dk[ell] <- (knotsPlusPlus[ell + k] - knotsPlusPlus[ell]) / k
        ints[2L, ell] <- dk[ell] * sum(BInt[2L, ell:(nKnots + k1 - 1L)])
    }
    
    ## Then substract the integrals from the first knot to 'xmin'
    for (ell in 1L:k1) {
        ints[1L, ell] <- dk[ell] * sum(BInt[1L, ell:(nKnots + k1 - 1L)])
    }

    ints <- ints[2L, ] - ints[1L, ]
    ## print(ints)

    b <- runif(ncol(B))
    
    if (FALSE && (k > 1L)) {
        b[1:k1] <- 0
        b[ncol(B) - 0L:(k - 2L)] <- 0
    }

    b <- b / (ints %*% b)   
    fp <- f  <- B %*% b
    
    ## df <- matrix(NA, nrow = length(knots), ncol = k)
    ## colnames(df) <- paste("ord", 0L:(k1))
  
    res <- list(x = x, f = f,
                xmin = xmin, xmax = xmax,
                order = k,
                ##leftDeriv = leftDeriv,
                ##rightDeriv = rightDeriv,
                knots0 = knots, 
                knots = knotsPlus,
                B = B, b = b, fp = fp)
                ## ints = ints,
                ## qp = resqp,
                ##R = R, r = r)

    class(res) <- "SplineDensity"
    res
    
    
}

