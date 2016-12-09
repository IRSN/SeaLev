\name{Brest.tide}
\Rdversion{1.1}
\alias{Brest.tide}
\docType{data}
\title{
  Density for the tidal level at Brest
}
\description{

  Gridded density for the tidal part of the sea level at Brest.

}

\usage{data(Brest.tide)}

\encoding{UTF-8} 

\format{

  The format is:
  List of 2
  \itemize{
    \item{$x}{ num [1:512] \code{100 101 102 102 103} ...}
    \item{$y}{ num [1:512] \code{0.00 5.84e-06 1.12e-05 1.66e-05 2.20e-05} ...}
  }
  
}

\details{

  This dataset gives an example of structure for further tidal
  densities. A desirable feature is that the values of \code{y} at end
  points should be zero (or numerically such) so that the first and last
  values of \code{x} indicate the Lowest possible value (Lowest
  Astronomical High Tide) and the Highest Astronomical Tide (HAT).

  The length of the vectors could have the form \eqn{2^p+1}{2^p+1} with
  \eqn{p}{p} integer for speed considerations in convolution. However in
  the present package the vector \code{y} will usually be interpolated
  to a longer one and an optimized choice of length will have no effect.
  
}

\source{
  
  The density was obtained by digitalizing the image on page 222 of
  B. Simon's book.

}

\references{

  B. Simon \emph{La \enc{marée}{maree} \enc{océanique}{oceanique}
    \enc{côtière}{cotiere}}.
  Institut \enc{océanographique}{oceanographique}, 2007.
  
}

\examples{
data(Brest.tide)
plot(Brest.tide, type = "l",
   lwd = 2, col = "SteelBlue3")
}
\keyword{datasets}
