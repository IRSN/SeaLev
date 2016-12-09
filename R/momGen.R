##*********************************************************************
## DO NOT ROXYGENIZE THIS!!! The .Rd files have been edited.
##*********************************************************************

##' Generating Function of the Moments or Cumulants.
##'
##' The generating function is computed by using a closed form
##' obtained by recursive integration by parts.
##' 
##' @title Generating Function of the Moments or Cumulants
##'
##' @param object An object of class \code{"SplineDensity"}.
##'
##' @param t Numeric vector of values at which the generating function
##' will be evaluated.
##' 
##' @param tmax When \code{t} is missing, a numerical vector
##' of values from \code{0.0} to \code{tmax} is used. 
##'
##' @param log Logical. If \code{TRUE} the reurned value is the
##' generating function of cumulants rather than moments.
##'
##' @return A numeric vector containing the values of the generating
##' function.
##'
##' @note The generating function of the cumulants often increases
##' rapidly with \eqn{t} and take very large values unless \eqn{t}
##' remains close to zero.
##' 
##' @details
##' The generating functions \eqn{M_x}{MX} and \eqn{K_X} for the moments 
##' are defined by \eqn{M_X(t) = \mathbb{E}[e^{tX}]}{M_X(t) = E[exp(tX)]} and
##' \eqn{K_X(t) = \log M_X(t)}{K_X(t) = log(M_X(t))}. The function \eqn{K_X}
##' is increasing and convex with \eqn{K_X(0) = 0}.
##'
##' See Wikipedia's page for
##' \href{https://en.wikipedia.org/wiki/Cumulant}{cumulants} and that for
##' \href{https://en.wikipedia.org/wiki/Moment-generating_function}{
##' moment generating function}.
##' 
##' @author Yves Deville
##'
##' @examples
##' opar <- par(mfrow = c(2, 1))
##' ck <- rSplineDens(order = 4, xmax = 10)
##' plot(ck, main = "spline density")
##' m <- momGen(ck, tmax = 2)
##' par(opar)
##'
##' 
momGen <- function(object, t = NULL, tmax = 0.4, log = TRUE) {

    check <- TRUE
    
    if (length(t) == 0) {
        t <- seq(from = 0, to = tmax, by = 0.02)
    }
    indtnn <- (t != 0.0)

    if (!any(indtnn)) {
        if (log) Mx <- rep(0.0, length(t))
        else Mx <- rep(1.0, length(t))
        return(list(x = t, y = Mx))
    }
    
    k <- object$order
    k1 <- k - 1L
    
    tau <- object$knots0
    tauPlus <- object$knots
    
    Btau <- splineDesign(knots = tauPlus,
                       x = tau,
                       ord = k,
                       derivs = rep(k1 , length(tau)),
                       outer.ok = TRUE)

    ftauDer <- Btau %*% object$b
    
    Mx <- rep(1, length(t))
    a <- diff(c(0, ftauDer))
    ## print(a)
    for (i in 1:length(t)) {
        if (t[i] != 0.0) {
            Mx[i] <- sum(a * exp(t[i] * tau)) / t[i]^k
        }
    }

    ## if m is odd, change sign
    if (k %% 2)  {
        Mx <- -Mx
    }

    ## correction for boundaries

    if (k >= 2L) {
        sign <- -1
        tauBound <- c(object$xmin, object$xmax)
        ## print(Mx)
        for (j in 0L:(k - 2L)) {
            fj <- splineDesign(knots = tauPlus,
                               x = tauBound,
                               ord = k,
                               derivs = rep(j, 2L),
                               outer.ok = TRUE) %*% object$b
            fj[2L] <- -fj[2L]
            ## cat(sprintf("derivative of order %d, %5.2f\n", j, fj))
            add <- sign * fj[1L] * exp(t[indtnn] * tauBound[1L]) / t[indtnn]^(j + 1)
            add <- add +  sign * fj[2L] * exp(t[indtnn] * tauBound[2L]) / t[indtnn]^(j + 1)
            sign <- -sign
            Mx[indtnn] <- Mx[indtnn] + add
        }
    }

    if (check && (length(t) > 2L)) {
        
        Mx.check <- rep(1, length(t))
        nx <- 256L
        x <- seq(from = object$xmin, to = object$xmax, length.out = nx)
        h <- x[2L] - x[1L]
        f <- predict(object, newdata = x)$y


        ## comput the first tow moments by trapezoidal rule. This
        ## could be done by a better method, but it is only for a
        ## wuick check.
        intd <- x * f
        Ex <- sum(intd[2:nx]) * h + (intd[1] + intd[nx]) * h / 2
        intd <- x * intd
        E2x <- sum(intd[2:nx]) * h + (intd[1] + intd[nx]) * h / 2
        Vx <- E2x - Ex^2

        cat(sprintf("Ex = %5.2f, Vx = %5.2f\n", Ex, Vx))

        ## plot(x, f, type = "l")
        
        for (i in 1:length(t)) {
            efi <- exp(t[i] * x) * f
            Mx.check[i] <- sum(efi[2:nx]) * h + (efi[1] + efi[nx]) * h / 2
        }

        plot(t, log(Mx), main = "generating function of cumulants and approx.",
             col = "orangered", type = "l", lwd = 2)
        abline(h = 0, v = 0, col = translude("gray", alpha = 0.6))
        Kapx <- t * (Ex + Vx * t / 2)
        lines(t, Kapx, col = "SteelBlue3", lty = "dashed", lwd = 2)
        abline(a = 0, b = Ex, col = "DarkOliveGreen3")
        lines(t, log(Mx.check), col = "orchid", lty = "dotted", lwd = 4)
        legend("topleft",
               legend = c("spline", "trap. rule", "1-st order", "2-nd order"),
               lwd = c(2, 2, 2, 2),
               lty = c("solid", "dotted", "solid", "dashed"),
               col = c("orangered", "orchid", "DarkOliveGreen3", "SteelBlue3"))
                   
    }

    if (log) Mx <- log(Mx)

    list(x = t, y = Mx)

}
