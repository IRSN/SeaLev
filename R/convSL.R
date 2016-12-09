##================================================================
## Author: Yves Deville
##
## Convolution for Sea Levels.
## 
##================================================================

convSLOld <- function(dens.x,
                      shift.x = 0,
                      threshold.y = NA,
                      distname.y = "exponential",
                      shift.y = ifelse(is.na(threshold.y), 0, threshold.y),
                      par.y = list(rate = 0.10),
                      covpar.y = NULL,
                      lambda = ifelse(is.na(threshold.y), 705.8, NA),
                      pct.conf = c(95, 70),
                      use.covlambda = "lambda" %in% colnames(covpar.y),
                      prob = NULL,
                      prob.max = ifelse(is.na(threshold.y), 1-1e-8, 1-1e-5),
                      pred.period = NULL,
                      N = 2048L,
                      plim.y = c(1e-12, 1-1e-12),
                      ## quad = c("trapmod", "trap"),
                      plot = TRUE,
                      show.x = TRUE,
                      show.asympt = TRUE,
                      alpha.below = 0.5,
                       trace = 0,
                      ...) {
    ## trace <- 0
    DEBUG <- TRUE
    
    ## Now the 'quad' arg is removed.
    quad = "trap"
    ## quad <- match.arg(quad)
    
    mc <- match.call()
    
    if (is.na(lambda)) stop("the rate 'lambda' must be given (in inverse years)")
    
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
    ## basic checks
    ##=================================================================
    
    nx <- length(dens.x$x)
    if (nx < 16) stop("at least 16 points must be given in 'dens.x'")
    
    ## checks for density
    if (is.unsorted(dens.x$x)) stop("unsorted x in dens.x")
    if (any(is.na(dens.x$y))) stop("NA not allowed in dens.x$y")
    if (any(dens.x$y < 0)) stop("<0 values not allowed in dens.x$y")
    
    if (any(dens.x$y[c(1, nx)]> 1e-6))
        warning("dens.x$y should be 0 at both end points")
    
    ##==================================================================
    ## Distribution names are as in 'fRenouv' function. Determination
    ## of the range of x and y
    ## ==================================================================
    
    special <- TRUE
    
    if (distname.y == "exponential") {
        funname.y <- "exp"
        parnames.y <- "rate"
    } else if (distname.y == "weibull") {
        funname.y <- "weibull"
        parnames.y <- c("shape", "scale")
    } else if (distname.y == "gpd") {
        funname.y <- "gpd"
        parnames.y <- c("scale", "shape")
    } else if (distname.y %in% c("log-normal", "lognormal")) {
        distname.y <- "log-normal"
        funname.y <- "lnorm"
        parnames.y <- c("meanlog", "sdlog")
    } else if (distname.y == "gamma") {
        funname.y <- "gamma"
        parnames.y <- c("shape", "scale")
    } else if (distname.y %in% c("MixExp2", "mixexp2")) {
        distname.y <- "mixexp2"
        funname.y <- "mixexp2"
        parnames.y <- c("prob1", "rate1", "delta")
    } else {
        
        parnames.y <- names(par.y)
        special <- FALSE
        
        warning("warning: distribution not in target list. Still EXPERIMENTAL")
        funname.y <- distname.y
        
    }
    
    ##====================================================================
    ## Wrapper functions
    ## 
    ## In this part, formals of functions are changed following the
    ## ideas in "fitdistr" of the MASS package. See therein.
    ##====================================================================
    
    ## reorder arguments to densfun and co
    dfun.y <- get(paste("d", funname.y, sep = ""), mode = "function")
    pfun.y <- get(paste("p", funname.y, sep = ""), mode = "function")
    qfun.y <- get(paste("q", funname.y, sep = ""), mode = "function")
    
    fms <- formals(dfun.y)
    args <- names(fms)
    m <- match(parnames.y, args)  
    if(any(is.na(m)))
        stop("'parnames.y' specifies names which are not arguments to 'densfun'")
    
    ## 'x' is maintened in pole position '1', then come the params in parnames,
    ## then the other if any (e.g. surrogate parameters)
    formals(dfun.y) <- c(fms[c(1L, m)], fms[-c(1L, m)])
    
    ## Caution: in function call, remember to use log = TRUE
    
    f.y <- function(parm, x) dfun.y(x, parm)
    
    ## Same thing for 'pfun'. Note that although the main arg of
    ## distributions functions is usually named "q", we force it to be
    ## "x" here, because 'f.y' and 'F.y' are used in the same manner
    ## in the log-likelihood!
    
    fms <- formals(pfun.y)
    args <- names(fms)
    m <- match(parnames.y, args)
    if(any(is.na(m)))
        stop("parnames.y specifies names which are not arguments to 'pfun'")
    
    formals(pfun.y) <- c(fms[c(1L, m)], fms[-c(1L, m)])
    
    F.y <- function(parm, x) pfun.y(x, parm)
    
    ## reorder arguments to densfun
    fms <- formals(qfun.y)
    args <- names(fms)
    m <- match(parnames.y, args)
    if(any(is.na(m)))
        stop("parnames.y specifies names which are not arguments to 'qfun'")
    
    formals(qfun.y) <- c(fms[c(1L, m)], fms[-c(1L, m)])
    
    q.y <- function(parm, p) qfun.y(p, parm)
    
    par.y <- unlist(par.y)
    p.y <- length(par.y)
    str <- paste("parm[", 1L:p.y, "]", collapse = ", ", sep = "")
    parnames.y <- names(par.y)
    
    str <- paste(paste(paste(parnames.y,
                             paste("parm[\"", parnames.y, "\"]", sep = ""),
                             sep = " = "),
                       collapse = ", "), sep = "")
    
    body(f.y) <- parse(text = paste("dfun.y(x,", str, ")"))     
    body(F.y) <- parse(text = paste("pfun.y(x,", str, ")"))
    body(q.y) <- parse(text = paste("qfun.y(p,", str, ")"))
    
    ##=========================================================================
    ## grid computations and evaluations
    ##
    ## WARNING: the convolution is computed without threshold.
    ##=========================================================================
    
    REDUCE <- FALSE
    
    ## grid characteristic for x
    x.min  <- dens.x$x[1L]
    x.max  <- dens.x$x[nx]
    xr <- x.max - x.min
    
    ## The same for y
    ylims <- q.y(par.y, p = plim.y)
    
    y.min <- ylims[1L]
    y.max <- ylims[2L]
    yr <- y.max - y.min
    
    if (trace) {
        cat("x.min and x.max : ", c(x.min, x.max),"\n")
        cat("y.min and y.max : ", c(y.min, y.max),"\n")
        cat("x.max - xmin :    ", xr, "\n")
        cat("y.max - ymin :    ", yr, "\n")
    }
    
    ## grid step
    h <- max(xr, yr) / N
    
    xg  <- seq(from = x.min, by = h, length.out = N)
    yg  <- seq(from = y.min, by = h, length.out = N) 
    
    z.min <- x.min + y.min + shift.x + shift.y
    zg <- seq(from = z.min, by = h, length.out = 2 * N - 1L)
    
    if (trace) {
        cat("z.min and z.max\n")
        if (REDUCE) print(c(z.min, zg[N]))
        else print(c(z.min, zg[2*N-1]))
    }
    
    ## Compute fX and FX.
    ## fX is obtained by interpolation
    fxg <- approx(x = dens.x$x, y = dens.x$y, xout = xg)$y
    fxg[is.na(fxg)] <- 0
    
    ## (re)normalize
    fxg <- fxg / sum(fxg) / h
    
    ## Fxg is only on used at return time
    Fxg <- cumsum(fxg) * h
    
    ## Compute fY and FY
    ## fyg <- f.y(par.y, x = yg)
    Fyg <- F.y(par.y, x = yg)
    
    if (quad == "trap") {
        fyg <- diff(c(Fyg, 1.0)) / h
        ## fyg <- f.y(par.y, x = yg)
    } else if (quad == "trapmod") {
        fyg <- diff(c(0.0, Fyg, 1.0), lag = 2) / 2 / h
    } else stop("bad value given for 'quad'")
    
    fzg <- convolve(fxg, rev(fyg), type = "o") * h
    condExp.x <- convolve(xg * fxg, rev(fyg), type ="o") * h / fzg
    condExp.x[1L] <- x.min
    
    Fzg <- cumsum(fzg) * h
    
    ## Change of 2010-09-15
    if (REDUCE) {
        zg <- zg[1L:N]
        fzg <- fzg[1L:N]
        Fzg <- Fzg[1L:N]
        condExp.x <- condExp.x[1L:N]
    }
    
    Tzg <- 1.0 / lambda / (1 - Fzg)
    
    ##==============================================================
    ## DELTA METHOD
    ## matrices for numerical derivation
    ## The column ip of Parmat contains the parameter value with
    ## a tiny modification of its ip component.
    ##==============================================================
    
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
        
        eps <- sqrt(.Machine$double.eps)
        ##eps <- 0.0001
        
        dparms <- abs(par.y) * eps
        dparms[dparms < eps] <- eps
        
        dFyg <- matrix(NA, nrow = N, ncol = p.y)
        dFzg <- matrix(NA, nrow = length(Fzg), ncol = p.y)
        dzzg <- matrix(NA, nrow = length(Fzg), ncol = p.y)
        
        colnames(dFyg) <- parnames.y
        colnames(dFzg) <- parnames.y
        colnames(dzzg) <- parnames.y
        
        for (ip in 1L:p.y) {
            par.plus <- par.y
            par.plus[ip] <- par.plus[ip] +  dparms[ip]
            par.minus <- par.y
            par.minus[ip] <- par.minus[ip] - dparms[ip]
            
            dFyg[ , ip] <- ( F.y(par.plus, x = yg) - F.y(par.minus, x = yg) ) / (2 * dparms[ip] )
            
            dProv <- convolve(fxg, rev(dFyg[ , ip]), type ="o")*h
            if (REDUCE) dFzg[ , ip] <- dProv[1:N]
            else dFzg[ , ip] <- dProv
            dzzg[ , ip] <- -dFzg[ , ip]/ fzg
        }
        
        sig.Fyg <- sqrt(apply(dFyg, 1, function(x){t(x) %*% covpar.yy %*% x }))
        sig.Fzg <- sqrt(apply(dFzg, 1, function(x){t(x) %*% covpar.yy %*% x }))
        
        if (use.covlambda) {
            cn.old <- colnames(dzzg)
            dzzg <- cbind((1 - Fzg) / lambda / fzg, dzzg)
            colnames(dzzg) <- c("lambda", cn.old)
            sig.zzg <- sqrt(apply(dzzg, 1, function(x){t(x) %*% covpar.y %*% x }))
        } else { 
            sig.zzg <- sqrt(apply(dzzg, 1, function(x){t(x) %*% covpar.yy %*% x }))
        }
        
    } else {
        
        dFyg <- matrix(NA, nrow = N, ncol = p.y)
        dFzg <- matrix(NA, nrow = length(Fzg), ncol = p.y)
        dzzg <- matrix(NA, nrow = length(Fzg), ncol = p.y)
        
        ## no uncertainty
        sig.Fyg <- rep(NA, length(Fzg))
        sig.Fzg <- rep(NA, length(Fzg))
        sig.zzg <- rep(NA, length(Fzg))
        
    }
    
    ##============================================================================
    ## Prepare results for RL plot
    ##============================================================================
    
    nc <- 3 + 2*length(pct.conf)
    cnames <-
        c("prob", "period", "quant",
          paste(rep(c("L", "U"), length(pct.conf)), rep(pct.conf, each = 2), sep = "."))
    
    ## Return levels (e.g. to be used in RLplot)
    ret.lev <- matrix(NA, nrow = length(prob), ncol = nc)
    rownames(ret.lev) <- prob
    colnames(ret.lev) <- cnames
    
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
            alpha.conf <- (100 - pct.conf[ipct])/100
            z.conf <- qnorm(1 - alpha.conf/2)
            ret.lev[ , 2*ipct + 2] <- zzg.app - z.conf * sig.zzg.app
            ret.lev[ , 2*ipct + 3] <- zzg.app + z.conf * sig.zzg.app
        }
        
    }
    
    ret.lev <- as.data.frame(ret.lev)
    
    ##============================================================================
    ## Prepare a pred matrix 
    ## Return levels (e.g. to be used in RLplot)
    ## Note that pred.prob is computed from pred.period
    ##
    ##============================================================================
    
    pred.prob <- 1 - 1 / lambda / pred.period
    ind <- (pred.prob>0) & (pred.prob < 1)
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
            alpha.conf <- (100 - pct.conf[ipct])/100
            z.conf <- qnorm(1 - alpha.conf/2)
            pred[ , 2*ipct + 2] <- zzg.app - z.conf * sig.zzg.app
            pred[ , 2*ipct + 3] <- zzg.app + z.conf * sig.zzg.app
        }
    }
    
    if (plot) {
        
        RSLplot(data = ret.lev,
                lambda = lambda,
                pct.conf = pct.conf,
                below = x.max + threshold.y,
                alpha.below = alpha.below,
                ...)
        
        abline(h = x.max, col = "gray")
        mtext(side = 4, at = x.max, text = "HAT", col = "gray")
        
        if (!is.na(threshold.y)) {
            abline(h = x.max + threshold.y, col = "red")
            mtext(side = 4,
                  at = x.max + threshold.y,
                  text = "HAT + thresh.",
                  col = "red")
        }
        
        if (show.x) {
            lines(x = -log( lambda * (1 - Fzg)),
                  y = condExp.x,
                  col = "purple")
        }
        
    }
    
    ## to be improved later 
    if (distname.y == "exponential") {
        Esp0 <- sum(exp(xg*par.y["rate"]) * fxg)*h
        Esp1 <- sum(xg*exp(xg*par.y["rate"]) * fxg)*h
        lmomExp <- c(log(Esp0), log(Esp1))
        mu.z <- threshold.y + log(Esp0) / par.y["rate"]
        sigma.z <- 1 / par.y["rate"]
        theopar.z <- c(mu.z, sigma.z)
        names(theopar.z) <- c("loc", "scale")
        
    } else {
        lmomExp <- c(NA, NA)
        theopar.z <- list()
    }
    
    ##================================================
    ## Plotting the asympotic behaviour
    ##================================================
    
    if (plot && show.asympt) {
        
        if (distname.y == "exponential") {    

            plot.asympt <- TRUE
            lRzg.asympt <- -(zg - shift.y)*par.y["rate"]  + log(Esp0)
            Fzg.asympt  <- 1.0 - exp(lRzg.asympt)
            T.asympt  <-  1.0 / lambda / (1.0 - Fzg.asympt)
            
        } else if (tolower(distname.y) == "gpd") {
            
            if (par.y["shape"] >= 0) {
                plot.asympt <- TRUE
                
                ## Approximatively exponential
                if (abs(par.y["shape"]) < 1e-4) {
                    if (trace) {
                        cat("GPD close to the exponential -> used in 'show.asympt'")
                    }
                    Esp0      <-  sum(exp(xg / par.y["scale"]) * fxg) * h
                    lRzg.asympt <- -(zg - shift.y) / par.y["scale"]  + log(Esp0)
                    Fzg.asympt  <- 1.0 - exp(lRzg.asympt)
                    T.asympt  <-  1.0 / lambda / (1.0 - Fzg.asympt)
                    
                } else {
                    
                    plot.asympt <- TRUE
                    delta  <- sum(xg*fxg)*h
                    ##cat("XXX delta", delta, "\n")
                    
                    Fzg.asympt  <-  pGPD(zg - delta,
                                         loc = shift.y,
                                         scale = par.y["scale"],
                                         shape = par.y["shape"])
                    
                    T.asympt  <-  1 / lambda / (1 - Fzg.asympt)
                    
                }
                
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
    
    
    ## n.max <- round(1.2*N)
    n.max <- 2 * N - 1L 
    
    res <- list(call = mc,
                threshold.y = threshold.y,
                shift.x,
                shift.y = shift.y,
                dens.x = list(x = xg, y = fxg),
                dist.x = list(x = xg, y = Fxg),
                dens.y = list(x = yg + shift.y, y = fyg),
                dist.y = list(x = yg +shift.y, y = Fyg),
                dens.z = list(x = zg[1L:n.max], y = fzg[1L:n.max]),
                dist.z = list(x = zg[1L:n.max], y = Fzg[1L:n.max]),
                z = zg[1L:n.max],
                Tz= Tzg[1L:n.max],
                logmomExp = lmomExp,
                condExp.x = condExp.x[1:n.max],
                theopar.z = theopar.z,
                dFy = dFyg,
                dFz = dFzg[1L:n.max, ],
                dzz = dzzg[1L:n.max, ],
                ret.lev = ret.lev,
                pred = as.data.frame(pred))
    
}

