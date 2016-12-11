\name{varnames.pdMat}
\title{
  
  GPtail variable names

}
\usage{
\method{varnames}{GPtail}(object, \dots)
\method{varnames}{GPtail}(object, \dots) <- value
}
\alias{varnames.GPtail}
\alias{varnames<-.GPtail}
\arguments{
  
  \item{object}{

    An object inheriting from class \code{GPtail} representing the
    convolution of a distribution with a spline density and a GPD.

  }
 
 \item{value}{

   A vector with the replacement values for the variable names
   associated with \code{object}. It must be a vector with length 3. The
   variables are: the variable with spline density distribution, the
   variable with GPD and the variable with the convolution distribution.


 }
 \item{\dots}{

   Not used in this method.

 } 
}
\description{
  
  These methods extract and set the coefficients associated with the
  model represented by  \code{object}. 

}
\value{
  Character vector of variable names.
}

\seealso{\code{\link{GPtail}}}

%% \examples{}

\keyword{models}
