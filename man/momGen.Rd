\name{momGen}
\alias{momGen}
\title{

  Generating Function of the Moments or Cumulants

}
\usage{

momGen(object, t = NULL, tmax = 0.4, log = TRUE)

}
\arguments{
  \item{object}{
    
    An object of class \code{"SplineDensity"}.
    
  }
  \item{t}{

    Numeric vector of values at which the generating function
    will be evaluated.

  }
  \item{tmax}{

    When \code{t} is missing, a numerical vector of values from
    \code{0.0} to \code{tmax} is used.

  }
  \item{log}{

    Logical. If \code{TRUE} the returned value is the generating function
    of cumulants rather than moments.

  }
  
}
\value{
  
A numeric vector containing the values of the generating function.

}
\description{

  Generating Function of the Moments or Cumulants.

}
\details{

  The generating function is computed by using a closed form obtained by
  recursive integration by parts.
  
  The generating functions \eqn{M_X} and \eqn{K_X} for the moments
  are defined by \eqn{M_X(t) = \mathrm{E}[e^{tX}]}{M_X(t) = E[exp(tX)]} and
  \eqn{K_X(t) = \log M_X(t)}{K_X(t) = log(M_X(t))}. The function \eqn{K_X}
  is increasing and convex with \eqn{K_X(0) = 0}.
  
  See Wikipedia's page for
  \href{https://en.wikipedia.org/wiki/Cumulant}{cumulants} and that for
  \href{https://en.wikipedia.org/wiki/Moment-generating_function}{
    moment generating function}.

}
\note{

  The generating function of the cumulants often increases
  rapidly with \eqn{t} and take very large values unless \eqn{t}
  remains close to zero.

}
\examples{
opar <- par(mfrow = c(2, 1))
ck <- rSplineDensity(order = 4, xmax = 10)
plot(ck, main = "spline density")
m <- momGen(ck, tmax = 2)
par(opar)
}
\author{
Yves Deville
}

