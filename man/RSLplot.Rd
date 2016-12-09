\name{RSLplot}
\alias{RSLplot}
\title{

  Return Sea Level plot

}
\encoding{UTF-8}
\description{
  
  Return Level Plot for Sea Level and convolution.

}
\usage{
   RSLplot(data,
           z = NULL,
           duration = NULL,
           lambda,
           xlab = "period", ylab = "level",
           mono  = TRUE,
           lty.quant = "solid",
           col.quant = ifelse(mono, "black", "SteelBlue"),
           pct.conf = c(95, 70),
           col.conf = NULL,
           filled.conf = FALSE,
           fill.conf = NULL,
           lty.conf = ifelse(rep(filled.conf, 2),
                             rep(NA, 2), c("dashed", "dotted")),
           rl.mark = 100,
           text.mark = rl.mark,
           col.mark = NULL,
           Tlim = NULL,
           problim = NULL,
           below = NA,
           alpha.below = 0.5,
           grid = TRUE,
           legend = TRUE,
           pch.points = c(24, 22, 23),
           col.points = c("black", "darkgray"),
           ...)
}
\arguments{
  
  \item{data}{

    Data frame containing the description of the return level curve and
    confidence limits. See \bold{Details}.

  }
  \item{z}{

    Empirical points to add to the plot.

  }
  \item{duration}{

    For empirical points.

  }
  \item{lambda}{

    Event rate (in events by year). This is used only when empirical
    points are added to the plot using the \code{x} argument.

  }
  \item{xlab}{

    Label for the x-axis.

  }
  \item{ylab}{

    Label for the y-axis.

  }
  \item{mono}{

    Monochrome plot.

  }
  \item{lty.quant}{

    Line type for the quantile (return level curve).

  }
  \item{col.quant}{

    Line color for the quantile (return level curve).

  }
  \item{pct.conf}{

    Confidence levels (in accordance with the content of the data.frame
    \code{data}.}  \item{col.conf}{Line colors for the confidence limits,
    in the same order as \code{pct.conf}.
    
  }
  \item{filled.conf}{

    Should confidence regions be filled?

  }
  \item{fill.conf}{

    Colors for the confidence regions if filled.

  }
  \item{lty.conf}{

    Line types for confidence regions limits. By using \code{NA}, the
    plot of the limits (lower and upper) will be cancelled which may be
    when \code{fill.conf} is \code{TRUE}.
    
  }
  \item{rl.mark}{

    Return level(s) to materialize with a vertical line.

  }
  \item{text.mark}{

    Label(s) to display near the \code{mark.rl} lines if any.

  }
  \item{col.mark}{

    Color(s) for the \code{marK.rl}.

  }
  \item{Tlim}{

    Limits for the x-axis given in time (return periods).

  }
  \item{problim}{

    Limits for the x-axis given in probability. Only one of the two
    arguments \code{Tlim} and \code{problim} can be given.

  }
  \item{below}{

    A level below which a semi-transparent rectangle will be drawn, with
    the meaning that the corresponding return levels should not be used.
    In the sea-level convolution context, this will normally be \code{HAT
      + threshold} because the distribution of the sea-level is conditional
    to a surge level exceeding the threshold.
    
  }
  \item{alpha.below}{

    A value of transparency to (partially) occult the
    region of the plot lying below the minimal level of validity for the
    convolution, see \bold{Details}. The value \code{0} means a fully
    transparent rectangle (no effect), and the value \code{1} means a
    fully opaque one. When using a device that do not support
    transparency, only the values \code{0} and \code{1} will be possible.
    
  }
  \item{grid}{

    Logical. If \code{TRUE}, a grid with dotted lines is added.
    
  }
  \item{legend}{
    
    Logical. If \code{TRUE}, a legend is added.
    
  }
  \item{pch.points}{

    Vector of symbols to be used for points given by \code{z} (if
    any). When \code{z} is a numeric vector, only the first element is
    used. When \code{z} is a list, the given symbols will be used in
    correspondence with the list elements.
       
  }
  \item{col.points}{

    Vector of colours used in the same way as \code{pch.points}.
       
  }
  \item{\dots}{

    Other arguments to pass to \code{plot}. These can be \code{main},
    \code{ylim}, etc.

  }
}
\details{
  
  This function is roughly a copy of the \code{RLplot} in the package
  \code{Renext} with a few supplementary formal arguments.
  
  Line types and colors for confidence regions are given as a vectors
  with elements corresponding to the percents given in
  \code{pct.conf}. These levels should be available in the data.frame
  \code{data}, which means that suitable columns should indicate the
  limits. For instance a \code{95} percent region can be displayed only
  when the two columns \code{L.95} and \code{U.95} exist and provide the
  relevant information.
  
}
\value{
 
  None.

}
%%\references{}
\author{

  Yves Deville

}
\seealso{

  The \code{RLplot} function of the \code{Renext} package.

}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
