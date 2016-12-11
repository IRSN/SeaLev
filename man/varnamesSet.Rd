\name{varnamesSet}
\alias{varnames<-}

\title{

  Assign names for variables

}

\usage{

varnames(object, \dots) <-  value

}

\alias{varnames<-}

\arguments{
  \item{object}{

    An object representing a structure with variables, e.g.  a
    multivariate model or a model involving several variables.

  }
  \item{\dots}{

    Some methods for this generic function may require
    additional arguments.

  }
  \item{value}{

    Value to be assigned to the variables associated with
    \code{object}. It must be a character vector with the
    same length as \code{varnames(object)}.

  }
}
\description{
  
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{GPtail}.

}
\value{
  
  Depends on the method function; see the appropriate documentation.

}

\seealso{\code{\link{varnames}}.}

%% \examples{ }
\keyword{models}
