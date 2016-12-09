##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************

##' Convolution of a spline density and a GP distribution.
##'
##' Given a density for a bounded r.v. \eqn{X}, the distribution of
##' \eqn{Z := X + Y} is computed where \eqn{Y} is a r.v. independent
##' of \eqn{X} with a Generalized Pareto distribution. The
##' r.v. \eqn{Y} will often be the excess in a POT model, in which
##' case the exceedance rate \eqn{lambda} should be used to get return
##' periods in years.
##' 
##' @title Convolution of a spline density and a GP distribution.
##'
##' @param x An object with S3 class \code{"SplineDensity"}.
##'
##' @param threshold.y Optional threshold.
##'
##' @param distname.y Name of the distribution in the Genezalized
##' Pareto family. For now, only \code{"GPD"} is allowed.
##'
##' @param par.y Numeric \emph{named} vector of GP parameters, see
##' examples. This is usually a vector of estimated parameters.
##'
##' @param covpar.y Covariance matrix to be used in a \emph{delta
##' method}. This is for the parameters of the \eqn{Y}-part, but can
##' contain a row and column for the event rate \code{lambda} (see
##' Details). The colnames and rownames must agree and must be equal
##' either to \code{names(par.y)} or to \code{c("lambda",
##' names(par.y))}.
##'
##' @param lambda Rate to be used in the computation of the return
##' levels. Should be given in \emph{events by year} since the return
##' levels are given on a yearly basis.
##'
##' @param shift.y
##' 
##' @param pct.conf Confidence levels in percent. Should be given in
##' decreasing order.
##'
##' @param use.covlambda Logical indicating if the uncertainty on the
##' event rate \code{lambda} should be taken into account in the delta
##' method or not.
##'
##' @param plot.which Integer choosing a plot, if any. The value
##' \code{0} plots nothing, \code{1} plot dentities and
##' \code{2} plot the survival..
##'
##' @param prob Probability for which the return levels are wanted in
##' the \code{ret.lev} table. A NULL value correspond to a default
##' vector of values.
##' 
##'
##' @param pred.period If not \code{NULL}, a vector giving periods at
##' which predictions (return levels and confidence limits) should be
##' computed and returned in the \code{ret.lev} data.frame.
##'
##' @param deriv Logical. If \code{TRUE}, the derivatives of the
##' survival of the sum \eqn{Z} will be computed.
##'
##' @param trace Integer level of verbosity.
##'
##' @param N Number of discretization points in the range of \code{x},
##' hence length of the vector used for \code{x} in the discrete
##' convolution.
##'
##' @param N.ex Number of (tail) values \eqn{z} for the r.v. q\eqn{Z}
##' where an exact computation of the survival \eqn{S_Z(z)} is
##' performed. The values are greater than the upper end-point of
##' \eqn{X}. These values are used in most results (density,
##' conditional expectation, ...) but not in the return levels
##' table where only rounded return periods are used after
##' an interpolation.
##'
##'
##' @return A list.
##' 
##' @author Yves Deville
##'
##' @examples
##' data(Brest.tide)
##' SD <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y,
##'                     order = 2, nKnots = 24)
##' par.y <- c("scale" = rgamma(1, shape = 2, scale = 30),
##'            "shape" = 0.2 * runif(1))
##' res <- GPtail(x = SD.Brest, par.y = par.y, lambda = 1)
##' 

GPtail <- function(x,
                   threshold.y = NA,
                   distname.y = "GPD",
                   par.y = c("scale" = 1, "shape" = 0.05),
                   covpar.y = NULL,
                   lambda = ifelse(is.na(threshold.y), 705.8, NA),
                   shift.y = ifelse(is.na(threshold.y), 0, threshold.y),
                   pct.conf = c(95, 70),
                   use.covlambda = "lambda" %in% colnames(covpar.y),
                   deriv = !is.null(covpar.y),
                   Tlim = c(1, 1e+5),
                   pred.period = NULL,
                   pred.prob = NULL,
                   trace = 1L,
                   N = 2 * 1024,
                   N.ex = 300,
                   plot.which = 0) {
    
    
    Ntot <- N + N.ex

    p.y <- 2L
    parnames.y <- c("scale", "shape")
    if (!setequal(parnames.y, names(par.y))) {
        stop("'par.y' must be a vector with names:", parnames.y) 
    }
    
    mc <- match.call()
    
    if (class(x) != "SplineDensity") {
        stop("'x' must be an object with class \"SplineDensity\".",
             " Use the 'SplineDensity' function to create such an object.")
    }
    
    if (distname.y != "GPD") {
        stop("For now, 'distname.y' can only be equal to \"GPD\"")
    }

    ##=================================================================
    ## default prob for quantile and confidence lims 
    ## should densify near 0 and 1
    ##=================================================================

    if (length(pred.period) && length(pred.prob)) {
        stop("only one of the formal 'pred.period' and 'pred.prob' can be given")
    }

    if (length(pred.prob)) {
        
        if (any(is.na(pred.prob))) stop("'pred.prob' values can not be NA") 
        if (any(pred.prob <= 0.0) || any(pred.prob >= 1.0)) {
            stop("'pred.prob' values must be > 0 and < 1")
        }
        pred.prob <- sort(pred.prob)
        npp <- length(pred.prob)
        pred.period <- 1.0 / (1.0 - pred.prob) / lambda
        if (pred.period[npp] > Tlim[2L]) {
            warning("'Tlim[2L]' increased to match prediction requirements")
            Tlim[2L] <- pred.period[npp]
        }
        
    } else if (length(pred.period)) {
        if (any(is.na(pred.period))) stop("'pred.period' values can not be NA")
        if (any(pred.prob <= 0.0)) {
            stop("'pred.period' values must be > 0")
        }
        pred.period <- sort(pred.period)
        pred.prob <- 1.0 - 1.0 / lambda / pred.period
    } else {
        rr <- ceiling(log(Tlim[2L], 10))
        if (rr > 4) { 
            pred.period <- (10^(rr - 4)) * c(0.1, 0.2, 0.5, 1:10, 20, 50, 100, 1000, 10000)
        } else {
            pred.period <- c(0.1, 0.2, 0.5, 1:10, 20, 50, 100, 1000, 10000)
        }
        pred.period <- pred.period[pred.period <= Tlim[2L]]
        pred.prob <- 1.0 - 1.0 / lambda / pred.period
    }

    if (trace) {
        cat("o Period for return levels (pred.period)\n")
        print(pred.period)
    }

    
    ##====================================================================
    ## manage the covariance matrix if needed
    ##====================================================================
    
    if (!is.null(covpar.y)) {
        
        covpar.y <- as.matrix(covpar.y)
        
        if ( (nrow(covpar.y) != ncol(covpar.y)) )
            stop("covariance matrix 'covpar.y' not square!")
        
        if ( is.null(colnames(covpar.y)) || is.null(colnames(covpar.y)))
            stop("matrix 'covpar.y' must have rownames and colnames")
        
        if (any(rownames(covpar.y) != rownames(covpar.y)))
            stop("matrix 'covpar.y' must have identical rownames and colnames")
        
        if ( !all(parnames.y %in% colnames(covpar.y)) ) {
            stop("colnames(covpar.y) does not contain all parnames")
        }
        
        parnames.all <- parnames.y

        if (use.covlambda) {
            if (!("lambda" %in% colnames(covpar.y))) {
                stop("when 'use.covlambda' is TRUE \"lambda\" must be",
                     " in the rownames and colnames of 'covpar.y'")
            }
            parnames.all <- c("lambda", parnames.all)
        }

        ## set row and column order
        covpar.y <- covpar.y[parnames.all, parnames.all, drop = FALSE]
        
    } else {
        use.covlambda <- FALSE
    }

    
    xmin <- x$xmin
    xmax <- x$xmax
    
    ymin <- 0
    if (par.y["shape"] < 0.0) ymax <- -par.y["scale"] / par.y["shape"]
    else ymax <- Inf

    ## XXX to be improved later.
    if (trace) {
        cat(sprintf("ymax = %6.1f\n", ymax))
        E.y <- par.y["scale"] / (1.0 - par.y["shape"])
        CV.y <- sqrt(1.0 / sqrt(1 - 2 * par.y["shape"]))
        sd.y <- E.y * CV.y
        cat(sprintf("E.y = %6.1f sd.y = %6.1f\n", E.y, sd.y))
    }
    
    Kx <- momGen(object = x, t = 1.0 / par.y["scale"], log = TRUE)$y
    muStar.x <- Kx * par.y["scale"]
    
    if (trace) {
        cat("o computing the expectation of the exponential tail: ")
        cat(sprintf("muStar.x = %6.3f\n", muStar.x))
    }
    
    ## chose grid points and eval density here.
    xg <- seq(from = xmin, to = xmax, length.out = N)
    h <- xg[2L] - xg[1L]
    fxg <- predict(x, newdata = xg)$y

    ##====================================================================
    ## chose the grid for Y. In fine, it will be the same as that
    ## for X. But we may want to check...
    ##====================================================================
    N.yg <- N
    yg <- seq(from = ymin, by = h, length = N.yg)
    if (par.y["shape"] < 0.0) {
        yg <- yg[yg < ymax]
        N.yg <- length(yg)
    }
    
    Fyg <- pGPD(q = yg , scale = par.y["scale"], shape = par.y["shape"])
    fyg <- diff(c(0.0, Fyg)) / h
    ## fyg <- f.y(par.y, x = yg)
    
    zmin <- xmin + ymin
    zg <- seq(from = zmin, by = h, length.out = N + N.yg - 1L)
    
    fzg <- convolve(fxg, rev(fyg), type = "o") * h
    Fzg <- cumsum(fzg) * h
    Szg <- 1.0 - Fzg
    
    ##====================================================================
    ## Compute the tail using the EXACT formula. This works for 'z'
    ## above 'xmax'
    ##====================================================================
    
    zmax <- xmax + qGPD(p = 1 / Tlim[2L] / lambda,
                        scale = par.y["scale"],
                        shape = par.y["shape"],
                        lower.tail = FALSE)

    if (trace) {
        cat(sprintf("o Using closed form smax = %8.0f\n", zmax))
    }
    
    zg.ex <- seq(from = xmax, to = zmax, length.out = N.ex)
    Szg.ex <- fzg.ex <- rep(0, N.ex)

    if (deriv) {
        dSzg.ex <- matrix(0, nrow = N.ex, ncol = 2L)
        colnames(dSzg.ex) <- c("scale", "shape")
    }
    
    k <- x$order
    k1 <- k - 1L
    tau <- x$knots0
    tauPlus <- x$knots
    
    ##=====================================================================
    ## Compute the coefficient 'a_i' and the border effects if needed.
    ## CAUTION a_j is stored in a[j + 1] for j = 0, 1, ...
    ##======================================================================
    xi.y <- par.y["shape"] 
    sigma.y <- par.y["scale"]
    
    if ((1.0 / xi.y) %in% 1L:k) {
        stop("the shape coefficient of the GPD must not be\n",
             " equal to an integer between 1 and the order")
    }
    
    a <- rep(NA, k)
    a[1L] <- (-sigma.y) / (1.0 - xi.y) 
    
    if (deriv) {
        da <- matrix(NA, nrow = k, ncol = 2L)
        colnames(da) <- c("scale", "shape")
        da[1L, "scale"] <- a[1L] / sigma.y 
        da[1L, "shape"] <- a[1L] / (1.0 - xi.y)
    } else {
        da <- NULL
    }

    ##===================================================================
    ## Unless k == 1, i.e. unless fx is a piecewise constant function,
    ## we need take into accound the k - 1 boundary corrections
    ##
    ## XXX TODO subset with 'indLeftPos' and/or 'infLeftPos' rather
    ## than add zeroes!
    ## ===================================================================
    
    if (k >= 2L) {
        
        tauBound <- c(xmin, xmax)
        
        ULeft <- 1.0 + xi.y * (zg.ex - xmin) / sigma.y
        SLeft <- ULeft^(-1.0 / xi.y)

        URight <- 1.0 + xi.y * (zg.ex - xmax) / sigma.y
        SRight <- URight^(-1.0 / xi.y)
        
        if (xi.y < 0.0) {
            indLeftPos <- (ULeft > 0.0)
            ULeft[!indLeftPos] <- SLeft[!indLeftPos] <- 0.0
            indRightPos <- (URight > 0.0)
            URight[!indRightPos] <- SRight[!indRightPos] <- 0.0
        } else {
            indLeftPos <- indRightPos  <- rep(TRUE, length(zg.ex))
        }
        
        ej <- 1.0 - xi.y
        
        for (j in 0L:(k - 2L)) {
            
            j1 <- j + 1L
            j2 <- j + 2L
            
            ## find derivatives at the two boundaries
            fj <- splineDesign(knots = tauPlus, x = tauBound,
                               ord = k,
                               derivs = rep(j, 2L),
                               outer.ok = TRUE) %*% x$b
            
            if (j == 0) atilde <- 1.0
            else atilde <- a[j]
            
            ## add term for 'fz' left (previous 'Sleft')
            add <- atilde * fj[1L] * SLeft
            ## add term for 'fz' right (previous 'Sright')
            add <- add  - atilde * fj[2L] * SRight
            fzg.ex <- fzg.ex - add
            
            ## add term for 'Sz' left 
            SLeft <- SLeft * ULeft
            add <- a[j1] * fj[1L] * SLeft
            ## add term for 'Sz', right
            SRight <- SRight * URight
            add <- add -  a[j1] * fj[2L] * SRight
            Szg.ex <- Szg.ex + add
            
            ## compute 'a' for next step
            a[j2] <- a[j1] * (-sigma.y) / (1 - j2  * xi.y)
            
            if (deriv) {
                
                ##====================================================
                ## Caution: SLeft and SRight are now for exponent
                ## -1.0 / xi.y + j + 1
                ## ===================================================
                
                BLeft <- ej * (1.0 - 1.0 / ULeft)
                BLeft[!indLeftPos] <- 0.0
                BRight <- ej * (1.0 - 1.0 / URight)
                BRight[!indRightPos] <- 0.0
                logULeft <- log(ULeft)
                logULeft[!indLeftPos] <- 0.0
                logURight <- log(URight)
                logURight[!indRightPos] <- 0.0
                
                ## derivative of 'Sz' w.r.t. 'xi', part 1
                add <- da[j1, "shape"] * fj[1L] * SLeft
                add <- add - da[j1, "shape"] * fj[2L] * SRight
                dSzg.ex[ , "shape"] <- dSzg.ex[ , "shape"] + add
                
                ## derivative of 'Sz' w.r.t. 'xi', part 2
                add <- a[j1] * fj[1L] * (logULeft - BLeft) * SLeft / xi.y / xi.y
                add <- add - a[j1] * fj[2L] * (logURight - BRight) * SRight / xi.y / xi.y
                dSzg.ex[ , "shape"] <- dSzg.ex[ , "shape"] + add
                
                ## derivative of 'Sz' w.r.t 'sigma', part 1
                add <- da[j1, "scale"] * fj[1L] * SLeft
                add <- add - da[j1, "scale"] * fj[2L] * SRight
                dSzg.ex[ , "scale"] <- dSzg.ex[ , "scale"] + add
                
                ## derivative of 'Sz' w.r.t 'sigma', part 2
                add <- a[j1] * fj[1L] * BLeft * SLeft / sigma.y / xi.y
                add <- add - a[j1] * fj[2L] * BRight * SRight / sigma.y / xi.y
                dSzg.ex[ , "scale"] <- dSzg.ex[ , "scale"] + add
                
                ## derivative of 'a' w.r.t 'xi' for next step
                rvec <- 1L:j2
                da[j2, "shape"] <- a[j2] * sum(rvec / (1 - rvec * xi.y))
                
                ## derivative of 'a' w.r.t 'sigma' for next step
                da[j2, "scale"] <- a[j2] * j2 / sigma.y

                ## for next step
                ej <- ej - xi.y
                
            }
            
        }
        
    }

    ## now take j = k - 1
    Btau <- splineDesign(knots = tauPlus,
                         x = tau,
                         ord = k,
                         derivs = rep(k1 , length(tau)),
                         outer.ok = TRUE)
    
    ftauDer <- Btau %*% x$b
    aDer <- a[k] * diff(c(0, ftauDer))
    
    if (k == 1L) atilde <- 1
    else atilde <- a[k1] 
    ## aPrimeDer <- atilde * diff(c(0, ftauDer))
    
    d <- diff(c(0, ftauDer))
    
    for (i in 1L:length(tau)) {
        
        Ui <- 1.0 + xi.y * (zg.ex - tau[i]) / sigma.y
        if (xi.y < 0.0) {
            indPos <- (Ui > 0.0)
            Ui <- Ui[indPos]
        } else {
            indPos <- rep(TRUE, N.ex)
        }
        
        Si <- Ui^(-1.0 / xi.y + k)
        
        Szg.ex[indPos] <- Szg.ex[indPos] + a[k] * d[i] * Si
        fzg.ex[indPos] <- fzg.ex[indPos] - atilde * d[i] * Si / Ui
        
        if (deriv) {
            
            Bi <- (1 - k * xi.y) * (1.0 - 1.0 / Ui)
            
            dSzg.ex[indPos, "shape"] <- dSzg.ex[indPos, "shape"] + da[k, "shape"] * d[i] * Si   
            dSzg.ex[indPos, "shape"] <- dSzg.ex[indPos, "shape"] + a[k] * d[i] *
                (log(Ui) - Bi) * Si / xi.y / xi.y

            dSzg.ex[indPos, "scale"] <- dSzg.ex[indPos, "scale"] + da[k, "scale"] * d[i] * Si   
            dSzg.ex[indPos, "scale"] <- dSzg.ex[indPos, "scale"] + a[k] * d[i] *
                Bi * Si / xi.y / sigma.y
            
        }
    }
    
    condExp.x.ex <-  (zg.ex + sigma.y / xi.y)  - Szg.ex / xi.y / fzg.ex

    if (deriv) {
        
        ## gradient of Syg
        dSyg <- matrix(NA, nrow = N.yg, ncol = 2L)
        colnames(dSyg) <- c("scale", "shape")
        Uyg <- 1.0 + xi.y * yg / sigma.y
        Syg <- 1.0 - Fyg
        
        Byg <- (1.0 - 1.0 / Uyg)
        dSyg[ , "scale"] <- Byg * Syg / xi.y / sigma.y 
        dSyg[ , "shape"] <- (log(Uyg) - Byg) * Syg / xi.y / xi.y        
        
        ## convolution
        dSzg <- matrix(NA, nrow = length(Fzg), ncol = 2L)
        colnames(dSzg) <- c("scale", "shape")
        
        dSzg[ , "scale"] <- convolve(fxg, rev(dSyg[ , "scale"]), type = "o") * h
        dSzg[ , "shape"] <- convolve(fxg, rev(dSyg[ , "shape"]), type = "o") * h
        
    } else {
        dSzg <- dSyg <- NULL
    }
    
    ##====================================================================
    ## Convolution
    ##====================================================================
   
    condExp.x <- convolve(xg * fxg, rev(fyg), type ="o") * h / fzg
    condExp.x[1L] <- xmin
    condExp.x[2L:4L] <- approx(x = xg[c(1L, 5L)],
                               y = condExp.x[c(1L, 5L)],
                               xout = xg[2L:4L])$y

    ##====================================================================
    ## Merge the two results: 'convolution' and 'exact'. Note that the
    ## discrete convolution approximation works for z < xmax since
    ## this corresponds to couples [x, y] below the main diagonal.
    ## ====================================================================
    
    ind <- (zg < xmax)
    zg <- c(zg[ind], zg.ex)
    fzg <- c(fzg[ind], fzg.ex)
    Szg <- c(Szg[ind], Szg.ex)
    Tzg <- 1.0 / lambda / Szg
    
    if (deriv) {
        dSzg <-  rbind(dSzg[ind, ], dSzg.ex)
        dzzg <- sweep(x = dSzg, MARGIN = 1L, STATS = fzg, FUN = "/")

        ## Note that for z -> xmin, the derivatives of z w.r.t the two
        ## parameters 'scale' and 'shape' tend to 0. However the
        ## derivative w.r.t 'lambda' tends to Inf(!), because the
        ## density of z is zero at z = xmin.
 
        for (j in 1L:2L) {
            dzzg[1L, j] <- 0
            dzzg[2L:4L, j] <- approx(x = zg[c(1L, 5L)],
                                     y = dzzg[c(1L, 5L), j],
                                     xout = zg[2L:4L])$y
        }
        if (use.covlambda) {
            dzzg <- cbind("lambda" = Szg / lambda / fzg, dzzg)
        } 
    } else {
        dzzg <- NULL
    }

    condExp.x <- c(condExp.x[ind], condExp.x.ex)
    
    ##====================================================================
    ## Compute standard deviations
    ##====================================================================
    
    if (!is.null(covpar.y)) {

        varFun <- function(x){
            t(x) %*% covpar.y[parnames.y, parnames.y] %*% x
        }
        varFun.all <- function(x){
            t(x) %*% covpar.y[parnames.all, parnames.all] %*% x
        }

        ## Fyg and Fzg do not depend on 'lambda'
        sig.Syg <- sqrt(apply(dSyg, 1, varFun))
        sig.Szg <- sqrt(apply(dSzg, 1, varFun))
        
        ## but the return level zzg depends on 'lambda'
        if (use.covlambda) {
            sig.zzg <- sqrt(apply(dzzg, 1, varFun.all))
        } else { 
            sig.zzg <- sqrt(apply(dzzg, 1, varFun))
        }

    } else {
        ## ignore uncertainty
        sig.Syg <- sig.Szg <- sig.zzg <- NULL
    }

    ##====================================================================
    ## SHIFT
    ##====================================================================

    yg <- yg + shift.y
    zg <- zg + shift.y
    
    ##====================================================================
    ## Compute return levels
    ##====================================================================

    prob <- c(## 0.0001,
        ## seq(from = 0.01, to = 0.09, by = 0.01),
        ## seq(from = 0.10, to = 0.80, by = 0.10),
        seq(from = 0.90, to = 0.99, by = 0.01),
        0.995, 0.996, 0.997, 0.998, 0.999,
        0.9995, 0.9996, 0.9997, 0.9998, 0.9999,
        0.99999, 0.999999, 0.9999999)
    prob.max <- 1.0 - 1.0 / lambda / Tlim[2L]
    prob <- prob[prob <= prob.max]
        
    Fzg <- 1 - Szg
    ind.nd <- !duplicated(Fzg)
    Fzg.nd <- Fzg[ind.nd]
    zg.nd <- zg[ind.nd]
    
    ret.lev <- matrix(NA, nrow = length(prob), ncol = 3L)
    rownames(ret.lev) <- NULL
    colnames(ret.lev) <- c("prob", "period", "quant")
    
    ret.lev[ , "prob"] <- prob
    ret.lev[ , "period"] <- 1.0 / lambda / (1.0 - prob)
    ret.lev[ , "quant"] <- approx(x = Fzg.nd, y = zg.nd, xout = prob)$y
    
    zzg.app <- approx(x = Fzg.nd, y = zg.nd, xout = prob)$y
    
    ret.conf <- matrix(NA, nrow = length(prob), ncol= 2L * length(pct.conf))
    colnames(ret.conf) <-
        paste(rep(c("L", "U"), length(pct.conf)), rep(pct.conf, each = 2L),
              sep = ".")
    if (!is.null(covpar.y)) {
        
        sig.zzg.app <- approx(x = zg, y = sig.zzg, xout = zzg.app)$y
        
        for (ipct in 1L:length(pct.conf)) {
            alpha.conf <- (100.0 - pct.conf[ipct]) / 100.0
            z.conf <- qnorm(1.0 - alpha.conf/2)
            ret.conf[ , 2L * ipct -1L] <- zzg.app - z.conf * sig.zzg.app
            ret.conf[ , 2L * ipct ] <- zzg.app + z.conf * sig.zzg.app
        }

    }
    
    ret.lev <- cbind(ret.lev, ret.conf)
    ret.lev <- as.data.frame(ret.lev)

    ##============================================================================
    ## Prepare a 'pred' matrix 
    ## Return levels (e.g. to be used in RLplot)
    ## Note that pred.prob is computed from pred.period
    ##
    ##============================================================================
    
    pred <- matrix(NA, nrow = length(pred.period), ncol = ncol(ret.lev))
    rownames(pred) <- pred.period
    colnames(pred) <- colnames(ret.lev)
    
    pred[ , "period"] <- pred.period
    pred[ , "prob"]   <- pred.prob
    pred[ , "quant"]  <- approx(x = Fzg.nd, y = zg.nd, xout = pred.prob)$y
    
    zzg.app <- approx(x = Fzg.nd, y = zg.nd, xout = pred.prob)$y
    
    if (!is.null(covpar.y)) {
        
        sig.zzg.app <- approx(x = zg, y = sig.zzg, xout = zzg.app)$y
        
        for (ipct in 1:length(pct.conf)) {
            alpha.conf <- (100 - pct.conf[ipct])/100
            z.conf <- qnorm(1 - alpha.conf/2)
            pred[ , 2L * ipct + 2L] <- zzg.app - z.conf * sig.zzg.app
            pred[ , 2L * ipct + 3L] <- zzg.app + z.conf * sig.zzg.app
        }
    }
    
    ##====================================================================
    ## plot if needed
    ##====================================================================

    if (plot.which) {
        
        if (plot.which == 1) {
            
            coll <- translude("gray", alpha = 0.9)
            
            opar <- par(mfrow = c(3, 1))
            par(mar = c(1, 4, 3, 3))
            
            plot(xg, fxg, type = "l", lwd = 2,
                 main = sprintf("sigma_Y = %6.2f, xi_Y = %6.2f",
                     par.y["scale"], par.y["shape"]),
                 col = "orangered", xlim = c(zmin, max(zg)))
            abline(h = 0, v = xmax, col = coll)
            par(mar = c(1, 4, 1, 3))
            plot(yg, fyg, type = "l", col = "SteelBlue3", lwd = 2)
            abline(h = 0, col = coll)
            par(mar = c(3, 4, 1, 3))
            plot(zg, fzg, xlim = c(zmin, max(zg)), type = "l",
                 col = "SpringGreen3", lwd = 2)
            abline(h = 0, v = xmax, col = coll)
            
            par(opar)
            
        } else if  (plot.which == 2) {
            coll <- translude("gray", alpha = 0.9)
            plot(zg, Szg, xlim = c(zmin, max(zg)), type = "l",
                 main = sprintf("sigma_Y = %6.2f, xi_Y = %6.2f",
                     par.y["scale"], par.y["shape"]),
                 col = "SpringGreen3", lwd = 2, ylim = c(0, 1))
            lines(zg.ex, Szg.ex, type = "l", lty = "dashed", lwd = 2)
            abline(h = 0, v = xmax, col = coll)
        } else if  (plot.which == 3) {
            RSLplot(data = x$ret.lev, lambda = x$lambda)
        }

    }
        
    res <- list(call = mc,
                lambda = lambda,
                threshold.y = threshold.y,
                shift.y = shift.y,
                SD.x = x,
                par.y = par.y,
                x = xg,
                dens.x = list(x = xg, y = fxg),
                y = yg,
                ## dist.x = list(x = xg, y = Fxg),
                dens.y = list(x = yg, y = fyg),
                dist.y = list(x = yg, y = Fyg),
                dens.z = list(x = zg, y = fzg),
                dist.z = list(x = zg, y = Fzg),
                deriv = deriv,
                a = a,
                da = da,
                aDer = aDer,
                z = zg,
                Tz = Tzg,
                condExp.x = list(x = zg, y = condExp.x),
                muStar.x = muStar.x,
                dSy = dSyg,
                dSz = dSzg,
                dzz = dzzg,
                ret.lev = ret.lev,
                pred = as.data.frame(pred))
    
    class(res) <- "GPtail"
    res
    
}


