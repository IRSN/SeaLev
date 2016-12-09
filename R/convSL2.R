##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************

##' Computes the distribution and return levels for sea levels from
##' the two parts 'tidal' and 'surge' using a convolution method.
##'
##' This function performs essentially the same thing as
##' \code{\link{convSL}}, but the numerical precision has been
##' improved by limiting the use of a discrete convolution to a
##' restricted range of values. The evaluation of the convolution tail
##' if obtained by using numerical quadratures for each of the tail
##' value. The improvements are for the case where the distribution of
##' \eqn{Y} is heavy-tailed and where large return periods are
##' considered, say > 1000 years.
##'
##' This function should replace the \code{\link{convSL}} function in
##' a future version of the package.
##' 
##' @title Convolution for Sea Levels: tide and surge
##' @param dens.x
##' @param shift.x 
##' @param threshold.y 
##' @param distname.y 
##' @param shift.y 
##' @param threshold.y
##' @param par.y 
##' @param covpar.y 
##' @param lambda 
##' @param pct.conf 
##' @param use.covlambda 
##' @param prob 
##' @param prob.max  
##' @param pred.period 
##' @param N
##' 
##' @param N.quad Number of (large) values \eqn{z} of the r.v. \eqn{Z}
##' for which a numerical quadrature will be used.
##' 
##' @param Tlim Numeric vector of length 2 giving the limits for the
##' return periods. Currently only the second element \code{Tlim[2]}
##' is used.
##'
##' @param deriv Logical. If TRUE, the derivative (w.r.t. parameters) of the
##' survival and the return level are computed and returned.
##'
##' @param plot
##' @param show.x 
##' @param show.asympt 
##' @param alpha.below 
##' @param trace 
##' @param ...
##'
##' @return A list similar to that returned by \code{\link{convSL}}.
##'
##' @note Slight changes have been made in the formal arguments list:
##' \code{N.quad}, \code{Tlim} and \code{deriv} are new formals while
##' the technical parameter \code{plim.y} was removed.
##' 
##' @author Yves Deville
##'
##' @examples
##' data(Brest.tide)
##' SD <- SplineDensity(x = Brest.tide$x, f = Brest.tide$y,
##'                     order = 2, nKnots = 24)
##' 
##' ## use a simulated sample to get a possible covariance matrix
##' ## (along with the estimated coefficients).
##' n <- 400
##' set.seed(123)
##' par.y0 <- c("scale" = rgamma(1, shape = 2, scale = 30),
##'             "shape" = 0.1 + 0.2 * runif(1))
##' par.y0 <- c("scale" = 58, "shape" = 0.28)
##' Y <- rGPD(n, scale = par.y0["scale"], shape = par.y0["shape"])
##' fit <- Renouv(x = Y, threshold = 0, effDuration = n,
##'               distname.y = "GPD", plot = FALSE)
##' par.y <- coef(fit)[c("scale", "shape")]
##' covpar.y <- vcov(fit)
##' lambda <- coef(fit)["lambda"]
##'
##' ## use the convSL function first
##' dev.set(2)
##' res1 <- convSL(dens.x = Brest.tide, distname.y = "GPD",
##'                covpar.y = covpar.y, lambda = lambda,
##'                N = 8 * 1024, par.y = par.y, plot = FALSE)
##' Tmax <- 100000
##' RSLplot(data = subset(res1$ret.lev, period <= Tmax),
##'         lambda = 1, main = "old")
##'
##' ## use the convSL2 function
##' if (length(dev.list()) < 2) dev.new()
##' dev.set(3)
##' res2 <- convSL2(dens.x = Brest.tide,
##'                 distname.y = "GPD",
##'                 covpar.y = covpar.y,
##'                 lambda = lambda,
##'                 N = 2 * 1024,
##'                 par.y = par.y,
##'                 plot = FALSE)
##'
##' RSLplot(data = subset(res2$ret.lev, period <= Tmax),
##'         lambda = 1,
##'         main = sprintf("new, scale = %5.2f shape = %5.2f",
##'                        par.y["scale"], par.y["shape"]))
##'
##' ## the GPtail function can be used for comparisons.
##' res3 <- GPtail(x = SD,
##'                par.y = par.y,
##'                lambda = lambda,
##'                covpar.y = covpar.y,
##'                N = 2 * 1014)
##'
##' cbind(subset(res1$ret.lev, subset = (period <= Tmax),
##'              select = c("period", "quant")),
##'       subset(res2$ret.lev, subset = (period <= Tmax),
##'              select = c("period", "quant")))
##' 
convSL <- function(dens.x,
                   shift.x = 0,
                   threshold.y = NA,
                   distname.y = "exponential",
                   shift.y = ifelse(is.na(threshold.y), 0, threshold.y),
                   par.y = c("rate" = 0.10),
                   covpar.y = NULL,
                   lambda = ifelse(is.na(threshold.y), 705.8, NA),
                   pct.conf = c(95, 70),
                   use.covlambda = "lambda" %in% colnames(covpar.y),
                   prob = NULL,
                   prob.max = ifelse(is.na(threshold.y), 1-1e-8, 1-1e-5),
                   pred.period = NULL,
                   N = 2048,
                   N.quad = NULL,
                   Tlim = c(1, 1e+5),
                   deriv = TRUE,
                   ## quad = c("trapmod", "trap"),
                   plot = TRUE,
                   show.x = TRUE,
                   show.asympt = TRUE,
                   alpha.below = 0.5,
                   trace = 0,
                   ...) {

    DEBUG <- TRUE
    plim.y = c(1e-12, 1-1e-12)
    period <- NULL ## to avoid WARNING at check due to in a later 'subset'.

    if (is.null(N.quad) || is.na(N.quad)) {
        N.quad <- round(8 * log(Tlim[2], base = 10))
    }
    
    mc <- match.call()
    
    if (is.na(lambda)) stop("the rate 'lambda' must be given (in inverse years)")

    if (!missing(par.y)) {
        par.y <- unlist(par.y)
    }

    
    ##=================================================================
    ## default prob for quantile and confidence lims 
    ## should densify near 0 and 1
    ##=================================================================
    
    if (is.null(prob)) {
        prob <- c(## 0.0001,
            ## seq(from = 0.01, to = 0.09, by = 0.01),
            ## seq(from = 0.10, to = 0.80, by = 0.10),
            seq(from = 0.90, to = 0.99, by = 0.01),
            0.995, 0.996, 0.997, 0.998, 0.999,
            0.9995, 0.9996, 0.9997, 0.9998, 0.9999,
            0.99999, 0.999999, 0.9999999)
        prob <- prob[prob <= prob.max]
    } else {
        if (any(is.na(prob))) stop("'prob' values can not be NA") 
        if ( any(prob <= 0.0) || any(prob >= 1.0) ) {
            stop("'prob' values must be > 0 and < 1")
        }
        prob <- sort(prob)
    }
    
    if (is.null(pred.period)) {
        rr <- 3
        pred.period <- (10^rr) * c(0.1, 0.2, 0.5, 1:10, 20, 50, 100, 1000, 10000)
    } else {
        if (any(is.na(pred.period))) stop("'pred.period' values can not be NA") 
        pred.period <- sort(pred.period)
    }
    
    ##=================================================================
    ## basic checks on dens.x
    ##=================================================================
    
    nx <- length(dens.x$x)
    if (nx < 16L) stop("at least 16 points must be given in 'dens.x'")
    
    ## checks for density
    if (is.unsorted(dens.x$x)) stop("unsorted x in dens.x")
    if (any(is.na(dens.x$y))) stop("NA not allowed in dens.x$y")
    if (any(dens.x$y < 0)) stop("< 0 values not allowed in dens.x$y")
    
    if (any(dens.x$y[c(1L, nx)] > 1e-6))
        warning("dens.x$y should be 0 at both end points")

    dens.x.unif <- (sd(diff(dens.x$x)) < 1e-6) 
    
    ##==================================================================
    ## Distribution names are as in 'Renext:::Renouv'
    ## function. Determination of the range of x and y
    ## ==================================================================

    myDist <- checkDist(distname.y = distname.y)
    funname.y <- myDist$funname.y
    distname.y <- myDist$distname.y
    
    if (myDist$special.y) {
        parnames.y <- myDist$parnames.y
        if (!setequal(parnames.y, names(par.y))) {
            stop("'par.y' must be a vector with names:", parnames.y) 
        }
        par.y <- par.y[parnames.y]
    } else {
        parnames.y <- names(par.y)
        special <- FALSE
        warning("warning: distribution not in target list. Still EXPERIMENTAL")
        funname.y <- distname.y
    }
    
    parnames.all <- c("lambda", parnames.y)
    ## parnb.y <- length(parnames.y)
    
    ##=========================================================================
    ## build probability functions and find the characteristics
    ## of the parameters. The two objects 'p.y' and 'fixed.y' are computed here
    ##=========================================================================
    
    myFuns <-  makeFuns(funname.y = funname.y,
                        parnames.y = parnames.y,
                        fixed.par.y = NULL,
                        trace = 0) 
    p.y <- myFuns$p.y
    
    ##=========================================================================
    ## grid computations and evaluations
    ##
    ## WARNING: the convolution is computed without threshold.
    ##=========================================================================
    
    REDUCE <- FALSE
    
    ## Grid characteristic for x
    xmin  <- dens.x$x[1L]
    xmax  <- dens.x$x[nx]
    xr <- xmax - xmin
    
    ## The same for y
    pmax = 1.0 - 1 / Tlim[2L] / lambda
    if (pmax < plim.y[2L]) plim.y[2L] <- pmax
    ylim <- myFuns$q.y(parm = par.y, p = plim.y)
    ymin <- ylim[1L]
    ymax <- ylim[2L]
    yr <- ymax - ymin
    
    ## zmax is larger than needed here
    zmax <- xmax + ymax

    if (trace) {
        cat("xmin and xmax : ", c(xmin, xmax), "\n")
        cat("ymin and ymax : ", c(ymin, ymax), "\n")
        cat("xmax - xmin :    ", xr, "\n")
        cat("ymax - ymin :    ", yr, "\n")
    }
    
    ##======================================================================
    ## grid design. The grid will be limited to 2 * xr to avoid
    ## numerical problems.
    ##======================================================================
    
    hg <- xr / N
    Nxg <- 2 * N
    xg  <- seq(from = xmin, by = hg, length.out = Nxg)

    ## XXX 
    Nyg <- 2 * N
    yg  <- seq(from = ymin, by = hg, length = Nyg)
    
    Nzg <- Nxg + Nyg - 1L
    
    zmin <- xmin + ymin + shift.x + shift.y
    zg <- seq(from = zmin, by = hg, length.out = Nzg)
    
    ##======================================================================
    ## Compute fX and FX.  fX is obtained by interpolation.
    ##=======================================================================
    
    fxg <- approx(x = dens.x$x, y = dens.x$y, xout = xg[1L:N])$y
    fxg[is.na(fxg)] <- 0.0
    fxg <- fxg / sum(fxg) / hg  ## (re)normalize
    Fxg <- cumsum(fxg) * hg
    fxg <- c(fxg, rep(0.0, N))
    Fxg <- c(Fxg, rep(1.0, N))
    
    ##======================================================================
    ## Compute fY and SY. Note that the density is unormalised because the
    ## upper end-point of Y is usually Inf
    ##======================================================================    

    Fyg <- myFuns$F.y(parm = par.y, x = yg)
    fyg <- diff(c(0.0, Fyg)) / hg
    Syg <- 1.0 - Fyg
    
    ##======================================================================
    ## Density fZ and conditional expectation E[ X | Z = z]. For small
    ## z the computation of E[ X | Z = z] is unreliable hence is
    ## replaced by linear interpolation.
    ##======================================================================
    
    fzg <- convolve(fxg, rev(fyg), type = "o") * hg
    Szg <- 1.0 - cumsum(fzg) * hg
    
    condExp.x <- convolve(xg * fxg, rev(fyg), type = "o") * hg / fzg
    condExp.x[1L] <- xmin
    condExp.x[2L:4L] <- approx(x = xg[c(1L, 5L)],
                               y = condExp.x[c(1L, 5L)],
                               xout = xg[2L:4L])$y
    
    ##=======================================================================
    ## reduce the number of values in z that can be used, because the
    ## upper half of them is unreliable
    ##=======================================================================
    
    Nzg <- Nxg
    fzg <- fzg[1L:Nzg]
    Szg <- Szg[1L:Nzg]
    condExp.x <- condExp.x[1L:Nzg]
    zg <- zg[1:Nzg]
    zgmax <- zg[Nzg]
 
    ##======================================================================
    ## Perform more computations, if required, now using one
    ## quadrature for each return period. Note than we can use the
    ## original data for density 'x' and 'f' since interpolation would
    ## not improve precision.
    ## ======================================================================    
    
    Tzgmax <- 1.0 / lambda / Szg[Nzg]    
        
    if (Tzgmax < Tlim[2L]) {
        
        zg.quad <- seq(from = zgmax + hg, to = zmax, length.out = N.quad)

        ## this was uncessfully experimented: perform quadratures at
        ## (approximated) quantiles corresponding to an uniform grid
        ## of probability.
        if (FALSE) {
            p1.quad <- myFuns$F.y(parm = par.y, x = zgmax + hg - xmax)
            p2.quad <- myFuns$F.y(parm = par.y, x = zmax - xmax)
            p.quad <- seq(from = p1.quad, to = p2.quad, length.out = N.quad)
            zg.quad <- xmax + myFuns$q.y(parm = par.y, p = p.quad)
        }
            
        fzg.quad <- Szg.quad <- condExp.x.quad <- rep(NA, N.quad)

        if (dens.x.unif) {
            xg.quad <- dens.x$x
            fxg.quad <- dens.x$y
            hg.quad <- dens.x$x[2L] - dens.x$x[1L]
            nx.quad <- nx
        } else {
            xg.quad <- xg[1:N]
            fxg.quad <- fxg[1:N]
            hg.quad <- hg
            nx.quad <- N
        }
        
        nx1.quad <-  nx.quad - 1L

        
        for (iq in 1L:N.quad) {

            zx.quad <- zg.quad[iq] - shift.x - shift.y - xg.quad
            ff <- fxg.quad * exp(myFuns$logf.y(parm = par.y, x = zx.quad))
            Sy.quad <- 1.0 - myFuns$F.y(parm = par.y, x = zx.quad)
            fS <- fxg.quad * Sy.quad
            ## XXX compute fx differently here?
            
            fzg.quad[iq] <- sum(ff[2L:nx1.quad]) * hg.quad +
                (ff[1L] + ff[nx.quad]) * hg.quad / 2

            Szg.quad[iq] <- sum(fS[2L:nx1.quad]) * hg.quad +
                (fS[1L] + fS[nx.quad]) * hg.quad / 2

            ff <- ff * xg.quad
            condExp.x.quad[iq] <- (sum(ff[2L:nx1.quad]) * hg.quad +
                                      (ff[1L] + ff[nx.quad]) * hg.quad / 2) /
                fzg.quad[iq]

        }

        ##===================================================================
        ## Merge the two results: 'convolution' and 'exact'.
        ##====================================================================

        comput <- c(rep("convol", length(zg)),
                    rep("quad", length(zg.quad)))
        zg <- c(zg, zg.quad)
        fzg <- c(fzg, fzg.quad)
        Szg <- c(Szg, Szg.quad)
        condExp.x <- c(condExp.x, condExp.x.quad)

    } else {
        N.quad <- 0L
        comput <- rep("convol", length(zg))
    }
    
    Tzg <- 1.0 / lambda / Szg
    
    ##===================================================================
    ## DELTA METHOD: do we need derivation w.r.t the rate 'lambda'?
    ##===================================================================
    
    if (!is.null(covpar.y)) {
        
        covpar.y <- as.matrix(covpar.y)
        
        if ( (nrow(covpar.y) != ncol(covpar.y)) )
            stop("covariance matrix 'covpar.y' not square!")
        
        if ( is.null(colnames(covpar.y)) || is.null(colnames(covpar.y)))
            stop("matrix 'covpar.y' must have rownames and colnames")
        
        if (any(rownames(covpar.y) != rownames(covpar.y)))
            stop("matrix 'covpar.y' must have identical rownames and colnames")
        
        if ( !all(parnames.y %in% colnames(covpar.y)) ) {
            stop("colnames(covpar.y) does not contain all par names")
        }
        
        covpar.yy <- covpar.y[parnames.y, parnames.y, drop = FALSE]

    }

    ##==================================================================
    ## DELTA METHOD: matrices for numerical derivation: 'par.plus' and
    ## 'par.minus' are p x p matrices with column j 'par.y +-' a small
    ## step in the direction of the j-th axis.
    ##===================================================================
    
    if (deriv) {

        eps <- sqrt(.Machine$double.eps)
        dSyg <- matrix(NA, nrow = Nyg, ncol = p.y)
        dSzg <- matrix(NA, nrow = Nxg + Nyg - 1L, ncol = p.y)
        
        colnames(dSyg) <- colnames(dSzg) <- parnames.y
  
        dparms <- abs(par.y) * eps
        dparms[dparms < eps] <- eps


        ## the 'if' is due to problems of lost names.
        if (p.y > 1L) {
            
            par.plus <- par.minus <- matrix(par.y, nrow = p.y, ncol = p.y)
            diag(par.plus) <- diag(par.plus) + dparms
            diag(par.minus) <- diag(par.minus) - dparms
            
            colnames(par.plus) <- colnames(par.minus) <- parnames.y
            rownames(par.plus) <- rownames(par.minus) <- parnames.y
            
            for (ip in 1L:p.y) {
                
                dSyg[ , ip] <- -(myFuns$F.y(parm = par.plus[ , ip], x = yg) -
                                     myFuns$F.y(parm = par.minus[ , ip], x = yg)) / (2 * dparms[ip] )
            
                dSzg[ , ip] <- convolve(fxg, rev(dSyg[ , ip]), type = "o") * hg            
                ## dzzg[ , ip] <- dSzg[ , ip] / fzg
                
            }
        } else if (p.y == 1L) {
            par.plus <- par.y + dparms
            par.minus <- par.y - dparms
            
            dSyg[ , 1L] <- -(myFuns$F.y(parm = par.plus, x = yg) -
                                     myFuns$F.y(parm = par.minus, x = yg)) / (2 * dparms)
            
            dSzg[ , 1L] <- convolve(fxg, rev(dSyg[ , 1L]), type = "o") * hg            
            ## dzzg[ , ip] <- dSzg[ , ip] / fzg
            
        }
            

        dSzg <- dSzg[1:Nxg, , drop = FALSE]
        
        ## compute the quadratures for the gradient we do not store
        ## a value for 'dSyg' here
        if (N.quad) {
            
            dSzg.quad <- matrix(NA, nrow = N.quad, ncol = p.y)
            
            for (iq in 1L:N.quad) {
                
                zx.quad <- zg.quad[iq] -shift.x - shift.y - xg.quad
                if (p.y > 1L) {
                    for (ip in 1L:p.y) {
                        
                        dS <- - (myFuns$F.y(par.plus[ , ip], x = zx.quad) -
                                     myFuns$F.y(par.minus[ , ip], x = zx.quad)) / (2 * dparms[ip])
                        fdS <- fxg.quad * dS
                        dSzg.quad[iq, ip] <- sum(fdS[2L:nx1.quad]) * hg.quad +
                            (fdS[1L] + fdS[nx.quad]) * hg.quad / 2
                    }
                } else if (p.y == 1L) {
                    dS <- - (myFuns$F.y(par.plus, x = zx.quad) -
                                 myFuns$F.y(par.minus, x = zx.quad)) / (2 * dparms)
                    fdS <- fxg.quad * dS
                    dSzg.quad[iq, 1L] <- sum(fdS[2L:nx1.quad]) * hg.quad +
                        (fdS[1L] + fdS[nx.quad]) * hg.quad / 2

                    
                }
                
            }
            
            ## Merge the two results: 'convolution' and 'exact'.
            dSzg <- rbind(dSzg, dSzg.quad)
            
        }

        ## compute the gradient of the return levels
        dzzg <- sweep(x = dSzg, MARGIN = 1, STATS = fzg, FUN = "/")
        
        if (!is.null(covpar.y)) {

            ## note that noe of the survivals in 'Syg' and 'Szg'
            ## depends on 'lambda' but the return level in 'zzg' does
            
            sig.Syg <- sqrt(apply(dSyg, 1, function(x){t(x) %*% covpar.yy %*% x }))
            sig.Szg <- sqrt(apply(dSzg, 1, function(x){t(x) %*% covpar.yy %*% x }))

            if (use.covlambda) {
                dzzg <- cbind("lambda" = Szg / lambda / fzg, dzzg)
                sig.zzg <- sqrt(apply(dzzg, 1, function(x){t(x) %*% covpar.y %*% x }))
            } else { 
                sig.zzg <- sqrt(apply(dzzg, 1, function(x){t(x) %*% covpar.yy %*% x }))
            }
            
        } else {
            ## no uncertainty
            sig.Syg <- sig.Szg <- sig.zzg <- NULL
        }
            
    } else {
        dSyg <- dSzg <- dzzg <- NULL
    }
    
    
    ##============================================================================
    ## Prepare results for RL plot
    ##============================================================================
    
    nc <- 3L + 2L * length(pct.conf)
    cnames <-
        c("prob", "period", "quant",
          paste(rep(c("L", "U"), length(pct.conf)), rep(pct.conf, each = 2), sep = "."))
    
    ## Return levels (e.g. to be used in RLplot)
    ret.lev <- matrix(NA, nrow = length(prob), ncol = nc)
    rownames(ret.lev) <- prob
    colnames(ret.lev) <- cnames
    
    Fzg <- 1.0 - Szg
    ind.nd <- !duplicated(Fzg)
    Fzg.nd <- Fzg[ind.nd]
    zg.nd <- zg[ind.nd]
    
    ret.lev[ , "prob"] <- prob
    ret.lev[ , "period"] <- 1 / lambda / (1 - prob)
    ret.lev[ , "quant"] <- approx(x = Fzg.nd, y = zg.nd, xout = prob)$y
    
    zzg.app <- approx(x = Fzg.nd, y = zg.nd, xout = prob)$y

    
    
    if (!is.null(covpar.y)) {
        
        sig.zzg.app <- approx(x = zg, y = sig.zzg, xout = zzg.app)$y
        
        for (ipct in 1:length(pct.conf)) {
            alpha.conf <- (100 - pct.conf[ipct]) / 100
            z.conf <- qnorm(1 - alpha.conf/2)
            ret.lev[ , 2 * ipct + 2] <- zzg.app - z.conf * sig.zzg.app
            ret.lev[ , 2 * ipct + 3] <- zzg.app + z.conf * sig.zzg.app
        }
        
    }
    
    ret.lev <- as.data.frame(ret.lev)
    ret.lev <- subset(ret.lev, period <= Tlim[2])
    
    ##============================================================================
    ## Prepare a pred matrix 
    ## Return levels (e.g. to be used in RLplot)
    ## Note that pred.prob is computed from pred.period
    ##============================================================================
    
    pred.prob <- 1 - 1 / lambda / pred.period
    ind <- (pred.prob > 0) & (pred.prob < 1)
    pred.period <- pred.period[ind]
    pred.prob <- pred.prob[ind]
    
    pred <- matrix(NA, nrow = length(pred.period), ncol = nc)
    rownames(pred) <- pred.period
    colnames(pred) <- cnames
    
    pred[ , "period"] <- pred.period
    pred[ , "prob"]   <- pred.prob
    pred[ , "quant"]  <- approx(x = Fzg.nd, y = zg.nd, xout = pred.prob)$y
    
    zzg.app <- approx(x = Fzg.nd, y = zg.nd, xout = pred.prob)$y
    
    if (!is.null(covpar.y)) {
        
        sig.zzg.app <- approx(x = zg, y = sig.zzg, xout = zzg.app)$y
        
        for (ipct in 1:length(pct.conf)) {
            alpha.conf <- (100.0 - pct.conf[ipct]) / 100.0
            z.conf <- qnorm(1 - alpha.conf/2)
            pred[ , 2 * ipct + 2L] <- zzg.app - z.conf * sig.zzg.app
            pred[ , 2 * ipct + 3L] <- zzg.app + z.conf * sig.zzg.app
        }
    }
    
    if (plot) {
        
        RSLplot(data = ret.lev,
                lambda = lambda,
                pct.conf = pct.conf,
                below = xmax + threshold.y,
                alpha.below = alpha.below,
                ...)
        
        abline(h = xmax, col = "gray")
        mtext(side = 4, at = xmax, text = "HAT", col = "gray")
        
        if (!is.na(threshold.y)) {
            abline(h = xmax + threshold.y, col = "red")
            mtext(side = 4,
                  at = xmax + threshold.y,
                  text = "HAT + thresh.",
                  col = "red")
        }
        
        if (show.x) {
            lines(x = -log(lambda * Szg),
                  y = condExp.x,
                  col = "purple")
        }
        
    }
    
    ## to be improved later 
    if (distname.y == "exponential") {
        Esp0 <- sum(exp(xg * par.y["rate"]) * fxg) * hg
        Esp1 <- sum(xg * exp(xg * par.y["rate"]) * fxg) * hg
        lmomExp <- c(log(Esp0), log(Esp1))
        mu.z <- threshold.y + log(Esp0) / par.y["rate"]
        sigma.z <- 1 / par.y["rate"]
        theopar.z <- c(mu.z, sigma.z)
        names(theopar.z) <- c("loc", "scale")
        
    } else {
        lmomExp <- c(NA, NA)
        theopar.z <- list()
    }
    
    ##============================================================
    ## Plotting the asympotic behaviour
    ##============================================================
    
    if (plot && show.asympt) {
        
        if (distname.y == "exponential") {    

            plot.asympt <- TRUE
            lRzg.asympt <- -(zg - shift.y) * par.y["rate"]  + log(Esp0)
            Fzg.asympt <- 1.0 - exp(lRzg.asympt)
            T.asympt <- 1.0 / lambda / (1.0 - Fzg.asympt)
            
        } else if (tolower(distname.y) == "gpd") {
            
            if (par.y["shape"] >= 0) {
                plot.asympt <- TRUE
                Esp0 <- sum(exp(xg / par.y["scale"]) * fxg) * hg
                
                plot.asympt <- TRUE
                delta  <- log(Esp0) * par.y["scale"]
                
                Fzg.asympt  <-  pGPD(zg - delta,
                                     loc = shift.y,
                                     scale = par.y["scale"],
                                     shape = par.y["shape"])
                
                T.asympt  <-  1.0 / lambda / (1.0 - Fzg.asympt)
                    
            } else {
                plot.asympt <- FALSE 
                warning("'show.asympt' not yet implemented for GPD with negative shape")
            }
        } else {
            plot.asympt <- FALSE 
            warning("'show.asympt' only for exponential or GPD distributions")
        }
        
        if (plot.asympt) {
            
            lines(x = log(T.asympt),
                  y = zg,
                  type = "l", lwd = 1, col = "SpringGreen3")
            
        }
    }
    
    res <- list(call = mc,
                threshold.y = threshold.y,
                shift.x,
                shift.y = shift.y,
                dens.x = list(x = xg, y = fxg),
                dist.x = list(x = xg, y = Fxg),
                dens.y = list(x = yg + shift.y, y = fyg),
                dist.y = list(x = yg + shift.y, y = 1.0 - Syg),
                dens.z = list(x = zg, y = fzg),
                dist.z = list(x = zg, y = 1.0 - Szg),
                z = zg,
                Tz = Tzg,
                zgmax = zgmax,
                Tzgmax = Tzgmax,
                condExp.x = condExp.x,
                logmomExp = lmomExp,
                theopar.z = theopar.z,
                dSy = dSyg,
                dSz = dSzg,
                dzz = dzzg,
                ret.lev = ret.lev,
                pred = as.data.frame(pred),
                comput = comput)
    
}
