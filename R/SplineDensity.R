##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************


##' Build a spline density from a provided grid density.
##'
##'
##' A spline approximation for a given density is found by a
##' constrained regression. First, a suitable basis of B-splines is
##' built. Then the coefficients are found in order to minimise the
##' distance to the provided density values with constraints arising
##' from boundary conditions and from the normalisation condition.
##' Boundary conditions can be given. By default, the values of the
##' density and of its first order derivative are taken equal to a
##' finite difference estimation from \code{x} and \code{f}. This
##' works correctly when the grid is fine enough, and when the
##' provided values correspond to those of a continuous function with
##' continuous derivative on the closed interval with end-points
##' \code{xmin} and \code{xmax}.
##' 
##' @title Build a spline density from a provided grid density
##'
##' @param x Numeric vector of values at which the density is provided.
##' 
##' @param f Numeric vector of density values corresponding to \code{x}.
##'
##' @param xmin Left (or lower) end-point of the distribution
##' 
##' @param xmax Right (or upper) end-point.
##'
##' @param leftDerivEst Integer vector giving the order of the
##' derivatives that will be estimated from the finite differences of
##' the points in \code{x} and \code{f}. These values will be used
##' even if derivatives are provided.
##'
##' @param rightDerivEst Similar to \code{leftDerivEst}. 
##'
##' @param leftDeriv Vector of known derivatives at left end-point (if
##' any). The given values are for order \eqn{0} to \eqn{m-2} in that
##' order. Unknown values are to be given as \code{NA}.
##' 
##' @param rightDeriv Similar to \code{rightDeriv} for the upper
##' end-point
##' 
##' @param knots A numeric vector of knots in ascending order.
##'
##' @param nKnots Number of knots to be used if \code{knots} is not
##' provided.
##'
##' @param order The spline order, e.g. \code{m = 4} for a cubic
##' spline.
##'
##' @param plot Logical. If \code{TRUE} a plot is provided.
##'
##' @param check Logical. When \code{TRUE}, some check of the
##' computations are carried over and results are printed.
##'
##' @return
##'
##' A list object that can be used for density computations. This object
##' is given an S3 class \code{"SplineDensity"}.
##'
##' @author Yves Deville
##'
##' @section Caution The spline is not warranted to be positive. 
##'
##' @seealso \code{\link{rSplineDensity}} to generate randomly drawn objects
##' (e.g. for tests).
##' 
##' @examples
##' 
##' data(Brest.tide)
##' SD <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y)
##' SD24 <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y, nKnots = 24)
##'
##' ## approximate a bounded GPD (negative shape) by a spline density
##' shape <- 2 + rexp(1)
##' x <- seq(from = 0, to = 1, length.out = 200)
##' f <- (1 - x )^(shape - 1) * shape
##' SDGP <- SplineDensity(x = x, f = f)
##' plot(SDGP)
##' 
SplineDensity <- function(x, f,
                          xmin = min(x),
                          xmax = max(x),
                          leftDerivEst = c(0L, 1L),
                          rightDerivEst = c(0L, 1L),
                          leftDeriv = NULL,
                          rightDeriv = NULL,
                          knots = NULL,
                          nKnots = 24L,
                          order = 4L,
                          plot = TRUE,
                          check = FALSE) {

    if (check) {
        checkInt <- checkConstr <- TRUE
    } else {
        checkInt <- checkConstr <- FALSE
    }
        
    nx <- length(x)
    k <- order
    k1 <- k - 1L

    if (!missing(xmin)) {
        if (xmin > min(x)) {
            stop("'xmin' must be <= min(x)")
        }
    }
    
    if (!missing(xmax)) {
        if (xmax < max(x)) {
            stop("'xmax' must be >= max(x)")
        }
    }
    xR <- xmax - xmin

    ##====================================================================
    ## Knots: no addition for boundary conditions at this stage
    ##====================================================================
    
    if (!missing(knots)) {
        if (knots[1] > xmin) knots <- c(xmin, knots)
        if (knots[length(knots)] < xmax) knots <- c(knots, xmax)
        if (!missing(nKnots))
            warning("'nKnots' ignored since 'knots' is given")
        nKnots <- length(knots)
    } else {
        knots <- seq(from = xmin, to = xmax, length.out = nKnots)
    }

    ##===================================================================
    ## Complete the sequence of knots. Note that w<e need 'k1' extra
    ## knots on both side. Yet we need one more on the right in order
    ## to compute the integrals later.
    ## ====================================================================

    knotsPlus <- knots
    eps <- 1e-3 * xR

    knotsPlus <- c(rep(knots[1] - eps, k1), knots, rep(knots[nKnots] + eps, k1))
    
    ##====================================================================
    ## Manage derivative information
    ##====================================================================
    
    ld <- length(leftDeriv)
    if (ld < k1) {
        if (!ld) leftDeriv <- rep(NA, k1)
        else leftDeriv <- c(leftDeriv, rep(NA, k1 - ld))
    }
    
    leftNm <- rightNm <- paste("ord ", 0:(k - 2L), ",", sep = "")
 
    ld <- length(rightDeriv)
    if (ld < k1) {
        if (!ld) rightDeriv <- rep(NA, k1)
        else rigthDeriv <- c(rightDeriv, rep(NA, k1 - ld))
    }

    ## Use finite difference if required for left derivatives
    if (0 %in% leftDerivEst) {
        leftDeriv[1] <- f[1L]
        leftNm[1] <- paste(leftNm[1], "from 'f'")
    }
    if (1 %in% leftDerivEst) {
        leftDeriv[2] <- (f[2L] - f[1L]) / (x[2L] - x[1L])
        leftNm[2] <- paste(leftNm[2], "from 'f'")
    }
    if (2 %in% leftDerivEst) {
        leftDeriv[3] <- (f[3L] - 2 * f[2L] + f[1L]) / (x[3L] - x[2L]) / (x[2L] - x[1L])
        leftNm[3] <- paste(leftNm[3], "from 'f'")
    }
    
    ## Use finite difference if required for right derivatives
    if (0 %in% rightDerivEst) {
        rightDeriv[1] <- f[nx]
        rightNm[1] <- paste(rightNm[1], "from 'f'")
    }
    if (1 %in% rightDerivEst) {
        rightDeriv[2] <- (f[nx] - f[nx- 1L]) / (x[nx] - x[nx - 1L])
        rightNm[2] <- paste(rightNm[2], "from 'f'")
    }
    if (2 %in% rightDerivEst) {
        rightDeriv[3] <- (f[nx] - 2 * f[nx - 1L] + f[nx - 2L]) / (x[nx] - x[nx - 1L]) /
            (x[nx - 1L] - x[nx - 2L])
        rightNm[3] <- paste(rightNm[3], "from 'f'")
    }

    names(leftDeriv) <- leftNm
    names(rightDeriv) <- rightNm
    
    cat("leftDeriv =", leftDeriv, "\n")
    cat("rightDeriv =", rightDeriv, "\n")
    
    ##====================================================================
    ## Make designs
    ##====================================================================

    ## 'n + k - 3' splines with indices '- k + 1' to 'n - 3'
    knotsPlusPlus <- c(knotsPlus, knotsPlus[length(knotsPlus)] + eps)
    BInt <- splineDesign(knots = knotsPlusPlus,
                         x = c(xmin, xmax),
                         ord = k + 1L,
                         outer.ok = TRUE)
    
    B <- splineDesign(knots = knotsPlus,
                      x = x,
                      ord = k,
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

    if (checkInt) {
        ## check by a numeric computation (trapezoidal rule)
        ints0 <- B[x >= xmin & x <= xmax, ]
        ni <- nrow(ints0)
        h <- diff(x[1L:2L])
        ints0 <- apply(ints0, 2, sum) * h - (ints0[1L, ] + ints0[ni, ]) * h / 2
        cat("\n")
        cat("o Check the integration of spline basis functions\n")
        cat("=================================================\n")
        print(cbind(coeff = round(dk, digits = 3),
                    numeric = round(ints0, digits = 3),
                    computed = round(ints, digits = 3),
              diff = round(ints - ints0, digits = 3)))
        
    }
    
    ##====================================================================
    ## build the matrix and vector of constaints
    ##====================================================================
    
    R <- matrix(rep(NA, 2 * k1 * ncol(B)), ncol = ncol(B))
    r <- rep(NA, 2 * k1)
    
    ## XXX FIXME bad length???
    ## rownames(R) <- names(r) <- c(paste("at min,", leftNm), paste("at max,", rightNm))

    if (k > 1L) {
        for (ii in 0L:(k - 2L)) {
            ind <- 1 + ii + c(0L, k1)
            
            R[ind, ] <- splineDesign(knots = knotsPlus,
                                     x = c(xmin, xmax),
                                     ord = k,
                                     derivs = rep(ii, 2),
                                     outer.ok = TRUE)
            r[ind] <- c(leftDeriv[ii + 1], rightDeriv[ii + 1])
        }
        
        ind <- !is.na(r)
        R <- R[ind, ]
        r <- r[ind]

    } 

    ## add normalisation constraint
    R <- rbind(R, ints)
    r <- c(r, 1.0)
        
    ##====================================================================
    ## solve for constraints
    ##====================================================================
   
    ## library(quadprog)

    resqp <- solve.QP(Dmat = t(B) %*% B,
                      dvec = t(B) %*% f,
                      Amat = t(R),
                      bvec = r,
                      meq = nrow(R))
    
    b <- resqp$solution

    if (checkConstr) {
        cat("\n")
        cat("o Check the constraints in quad. prog\n")
        cat("=====================================\n")
        print(cbind("Rb" = R %*% b, "r" = r, "diff" = unname(R %*% b - r)))
    }
        
    fp <- B %*% b

    if (any(fp < 0)) {
        warning("negative elements in 'fp' with min\n", min(fp[fp  < 0]))
        fp[fp < 0] <- 0
    }

    Bd <- splineDesign(knots = knotsPlus,
                       x = knotsPlus,
                       ord = k,
                       derivs = rep(k1 , length(knotsPlus)),
                       outer.ok = TRUE)

    res <- list(x = x, f = f,
                xmin = xmin, xmax = xmax,
                order = k,
                leftDeriv = leftDeriv,
                rightDeriv = rightDeriv,
                knots0 = knots, 
                knots = knotsPlus,
                B = B, b = b, fp = fp,
                ints = ints,
                qp = resqp,
                R = R, r = r)

    class(res) <- "SplineDensity"
    res
    
}
