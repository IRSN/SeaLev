\name{convSL}
\Rdversion{1.1}
\alias{convSL}
\title{
  
  Convolution for Sea Levels tide and surge

}
\encoding{UTF-8}
\description{
  
  Computes the distribution and return levels for sea levels from the
  two parts \sQuote{tidal} and \sQuote{surge} using a convolution
  method.

}
\usage{
   convSL(dens.x,
          shift.x = 0,
          threshold.y = NA,
          distname.y = "exponential",
          shift.y = ifelse(is.na(threshold.y), 0, threshold.y),
          par.y = c("rate" = 0.1),
          covpar.y = NULL,
          lambda = ifelse(is.na(threshold.y), 705.8, NA),
          pct.conf = c(95, 70),
          use.covlambda = "lambda" \%in\% colnames(covpar.y),
          prob = NULL,
          prob.max = ifelse(is.na(threshold.y), 1-1e-8, 1-1e-5),
          pred.period = NULL,
          N = 2048L,
          N.quad = NULL,
          Tlim = c(1, 1e+5),
          deriv = TRUE,
          plot = TRUE,
          show.x = TRUE,
          show.asympt = TRUE,
          alpha.below = 0.5,
          trace = 0,
          ...)
}
\arguments{
  
  \item{dens.x}{

    List containing the density for the tidal level.  It must contain
    two vectors \code{$x} and \code{$y} giving points on the density
    curve. The \code{x} vector must be sorted in ascending order, and
    the \code{y} vector must contain non negative values with end values
    equal to 0 or nearly such. The first and last values in \code{x}
    give the lowest and the highest astronomical tides up to the shift
    \code{shift.x}.

  } 
  \item{shift.x}{

    Numeric constant to add to the total level e.g.  a reference level
    if the density in \code{dens.x} has to be shifted. This can be
    useful to change the reference elevation.

  }
  \item{threshold.y}{

    Threshold used in the surge POT fit. The surge excesses over this
    threshold are assumed to follow a known distribution.

  }
  \item{distname.y}{
    
    The distribution name for the surge excesses over the threshold. The
    supported distribution names are those of the \code{fRenouv} in the
    \code{Renext} package. The distribution values are the non-negative
    reals. Most distributions will have one or two parameters: for a gpd
    distribution the, location parameter must be \code{0} since it
    applies to excesses.

  }
  \item{shift.y}{

    Numeric constant having the same dimension as the variable \code{Y},
    and used to shift the distribution. When a POT distribution is used
    this is automatically set to \code{threshold.y} and thus the
    distribution is assumed to be for the excesses over the
    threshold. In other cases, it should be normally zero.

  }
  \item{par.y}{

    Named list or vector containing the values of the
    parameters. Parameters with default values in the corresponding
    densities and distribution functions can be omitted and will be such
    if plugged in from a POT output. For instance the location parameter
    of the gpd will not be given and will take its default value
    \code{0}.

  }
  \item{covpar.y}{

    A covariance matrix for the parameters of the y-part. It can also
    contain a row and column for the event rate \code{lambda} (see
    Details). The colnames and rownames must agree and must be equal
    either to \code{names(par.y)} or to \code{c("lambda",
    names(par.y))}.

  }
  \item{lambda}{
    
    The event rate in the surge POT. Should be given in \emph{events by
    year} since the return levels are given on a yearly basis.

  }
  \item{pct.conf}{

    Confidence levels in percent. Should be given in decreasing order.

  }
  \item{use.covlambda}{

    Logical indicating if the uncertainty on the event rate
    \code{lambda} should be taken into account in the delta method or
    not.

  }
  \item{prob}{

    Probability for which the return levels are wanted in the
    \code{ret.lev} table. A \code{NULL} value correspond to a default
    vector of values.

  }
  \item{prob.max}{

    Maximal probability for the return level table/plot.

  }
  \item{pred.period}{

    If not \code{NULL}, a vector giving periods at which predictions
    (return levels and confidence limits) should be computed and
    returned in the \code{pred} data.frame. The results are returned in
    a \code{pred} data.frame.

  }
  \item{N}{

    Number of points in the convolution grid. Should be a power of two
    for a best computation speed.

  }
  \item{N.quad}{

    Number of quadrature points, corresponding to high return periods.

  }
  \item{Tlim}{

    Limits for the return periods.
    
  }
  \item{deriv}{

    Logical. If \code{TRUE} numerical derivatives are used to compute
    approximate confidence intervals by the delta method.
    
  }
  %% \item{quad}{
  %%   Quadrature rule underlying the convolution. Can
  %%   be \code{"trapmod"} or \code{"trap"}. The differences between
  %%   the two methods seem small in practice.
  %% }
  \item{plot}{
    
    Logical. Should a return level plot be drawn?
    
  }
  \item{show.x}{
    
    Logical indicating if a "tidal return level" curve should
    be added to the return level plot. See \bold{Details} below.
    
  }
  \item{show.asympt}{
    
    Logical. When \code{TRUE}, a curve showing the asymptotic behavior
    is added to the return level plot. This is available only for the
    exponential or the GPD with positive shape. In the first case, the
    added curve should be close to the curve obtained by convolution. In
    the second case, the added curve should be broadly parallel to the
    convolution curve. The added curve indicates the shape of the curve
    for \emph{very large} return periods.

  }
  \item{alpha.below}{

    A value of transparency to (partially) occult the region of the plot
    lying below the minimal level of validity for the convolution, see
    \bold{Details}. The value \code{0} means a fully transparent
    rectangle (no effect), and the value \code{1} means a fully opaque
    one. When using a device that do not support transparency, only the
    values \code{0} and \code{1} will be possible.

  }
  \item{trace}{

    Integer giving a level of verbosity. The value \code{0} leads to
    printing nothing.

  }
  \item{\dots}{

    Further arguments to be passed to \code{RSLplot} and then to
    \code{plot}. Most commonly used are \code{Tlim} \code{problim} from
    \code{RSLplot}, \code{main}, \code{ylim} from \code{plot}.  The
    arguments \code{z} an \code{duration} can also be used to add
    empirical points to the plot, see the vignette shipped with this
    package for examples.

  }
}
\details{
  
  The function computes the density and distribution functions for a
  surge \eqn{Z = X +Y}{Z = X +Y} where \eqn{X}{X} is a tidal level and
  \eqn{Y} is a random surge level with given (conditional) distribution
  over the threshold. The tidal level \eqn{X} has a distribution given
  in \code{dens.x}. The two parts \eqn{X} and \eqn{Y} are assumed to be
  independent hence the density of \eqn{Z}{Z} can be found using a
  numerical convolution.
  
  Since the distribution for the surge \eqn{Y}{Y} is given only for
  \eqn{Y > u }{Y > u} where \eqn{u}{u} is the threshold, the
  distribution of \eqn{Z}{Z} is known only for \eqn{Z > \textrm{HAT} +
  u}{Z > HAT + u}. Therefore \bold{only the corresponding part of the
  return level curve must be considered}. The left-side part can be
  considered as an extrapolation for small surge levels which is
  unwarranted.
 
  An approximate inference is derived using the "delta method" and the
  covariance matrix given in \code{covpar.y} (if any).  The return level
  \eqn{z(T)}{z(T)} corresponding to \deqn{z(T) = q_Z(p) \qquad p = 1 -
  \frac{1}{\lambda T}}{ z(T) = q.Z(p), p = 1 - 1 / lambda / T} where
  \eqn{q_Z(p)}{q.Z(p)} is the quantile function of \eqn{Z}{Z}.  The
  quantile function is deduced from the distribution function and the
  event rate \eqn{\lambda}{lambda} is replaced by its estimation
  provided in the \code{lambda} formal.
  
  Depending on the value of \code{use.covlambda} the impact on the
  estimation of \eqn{\lambda}{lambda} on the return level will be taken
  into account or not.
 
  When \code{plot} is \code{TRUE}, a return level plot is shown using a
  log scale for the return periods. The return level curve will have an
  asymptote for large return periods when the distribution of \code{Y}
  is exponential.

  When \code{plot} and \code{show.x} are \code{TRUE}, the conditional
  expectation curve is shown. This curve shows the conditional expected
  tidal level \deqn{\textrm{E}[X \vert Z=z]}{E[X|Z=z]} with
  \eqn{z=z(T)}{z=z(T)} taken equal to the return level associated to the
  period \eqn{T}{T}. This curve might not be visible and it can be
  necessary to adjust the \code{ylim} argument. The curve points out the fact
  that return levels for sea levels are \emph{not} conditional to the
  tidal level, and therefore that the sea levels corresponding to a
  given tide are much greater than suggested by the return level curve.

  It can be shown that when the distribution of surge is exponential
  the conditional expectation is constant for large return periods
  \eqn{T}{T} and then takes a value < HAT. When the distribution of surge is GPD
  with a positive shape,
  the conditional expectation is increasing but with very slow rate.
  
}
 
\value{
  
  A list with elements
  \item{dens.x, dens.y, dens.z}{

    Density of the tidal level, the surge and the total sea level.

  }
  \item{dist.x, dist.y, dist.z}{

    Distribution functions as for the densities

  }
  \item{z, T}{

    Numeric vectors of return levels and return periods. These
    are all the computed values and the values of \code{T} do not
    take only 'pretty' values. 

  }
  \item{dSy, dSz, dzz}{

    Matrices with columns giving the derivatives of \code{Sy}, \code{Sz}
    or \code{z} with respect to the parameters taken from
    \code{par.y}. For \code{dz} an extra column is placed in position
    \code{1} containing the derivative with respect to the event rate if
    \code{use.covlambda} is \code{TRUE}.

  } 
  \item{ret.lev}{

    Data.frame containing the return levels for some given or default
    probabilities. It contains columns for the probability, the return
    period, the return level (or quantile) as well as confidence limits
    for the return level. When \code{use.covlambda} is \code{TRUE} the
    confidence limits should be associated with the return periods (and
    not the probabilities), at least for relatively small probabilities.
    
  }
  \item{pred}{

    Data.frame containing the predicted return levels for some given or
    default return levels.

  }
  \item{condExp.x}{

    The conditional expectation vector, in correspondence with the grid
    for the density of \code{x}.

  }
  \item{logmomExp}{

    When the distribution of \eqn{Y}{Y} is exponential, the logarithm of
    the exponential moment
    \eqn{\textrm{E}[e^{X/\sigma_Y}]}{E[exp(X/sigma.y)]} where
    \eqn{\sigma_Y}{sigma.y} is the scale parameter of the exponential
    distribution. See the vignette shipped with this package for more
    information.

  }
  \item{logmomExp}{
    Character vector 
  }
  
}

\references{
  
  D. T. Pugh and J. M. Vassie (1978)
  "Extreme sea-levels from tide and surge probability"
  \emph{Proceeding of the 16th Coastal Engineering Conference}, Hamburg.
  
  B. Simon \emph{La \enc{marée}{maree} \enc{océanique}{oceanique}
    \enc{côtière}{cotiere}}.
  Institut \enc{océanographique}{oceanographique}, 2007.

}
\author{

  Yves Deville

}
\section{Warnings}{

  The densities and conditional expectations are given on for a wide
  grid and for very small probability of exceedance
  \eqn{1-F_Z(z)}{1-F.Z(z)}, say less than \code{1e-4}. The numerical
  precision weakens for very small values of \eqn{F_Z(z)}{F.Z(z)}. In
  future versions, the vectors will be truncated at a suitable value. A
  rule for a robust determination (with respect to the distribution of
  surges) still needs some investigations. Yet a minimal probability of
  exceedance about \code{1e-4} seems a good indication.
  
  Some distribution of the surge can have unbounded density,
   e.g. Weibull or gamma. This may affect the precision of the discrete
   approximation used in the numerical convolution.
  
}

\note{
  
  When a non-zero value is used for \code{shift.x}, the returned density
  for \eqn{X}{X} is left unchanged, i.e. is not shifted. Only the
  distribution of \eqn{Z} is shifted. On the contrary when an non-zero
  \code{shift.y} is given, the density of \eqn{Y}{Y} is shifted. Thus if
  a POT model is given with threshold \eqn{u}{u} the distribution
  specified on input is for the excesses \eqn{Y-u}{Y-u}, but the
  \code{dens.y} and \code{dist.y} objects of the returned list are for
  the \eqn{Y} levels.

}

\seealso{
  
  \code{\link{RSLplot}} for return level plot. The previous version
  \code{\link{convSLOld}}.
  
}

\examples{

## This example NEEDS Renext 
data(Brest); data(Brest.tide)

## POT
distname.y <- "exponential"
fit.exp <-
  Renouv(x = Brest$OTdata$Surge,
         effDuration = as.numeric(Brest$OTinfo$effDuration),
         threshold = 50, distname.y = distname.y,
         main = "exponential POT for surge")

## plug results into convSL
res.exp <-
   convSL(dens.x = Brest.tide,
          threshold.y = 50, distname.y = distname.y,
          lambda = fit.exp$estimate["lambda"],
          par.y = fit.exp$estimate["rate"],
          covpar.y = fit.exp$cov,
          use.covlambda = TRUE, main = "exponential")

## some results 
round(res.exp$pred, digits = 0)
plot(res.exp$dens.z, type ="l",
     main = "density of Z", xlab = "m", ylab = "") 

## POT with gpd with two
distname.y <- "gpd"

fit.gpd <-
   Renouv(x = Brest$OTdata$Surge,
          effDuration = as.numeric(Brest$OTinfo$effDuration),
          threshold = 50, distname.y = distname.y,
          main = "gpd POT for surge")
res.gpd1 <-
   convSL(dens.x = Brest.tide,
          threshold.y = 50, distname.y = distname.y,
          lambda = fit.gpd$estimate["lambda"],
          par.y = fit.gpd$est.y,
          covpar.y = fit.gpd$cov,
          use.covlambda = TRUE,
          main = "gpd with uncertainty on \"lambda\"")

## ignore the uncertainty on 'lambda'
res.gpd2 <-
    convSL(dens.x = Brest.tide,
           threshold.y = 50, distname.y = distname.y,
           lambda = fit.gpd$estimate["lambda"],
           par.y = fit.gpd$est.y,
           covpar.y = fit.gpd$cov,
           use.covlambda = FALSE,
           main = "gpd ignoring uncertainty on \"lambda\"")

}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
