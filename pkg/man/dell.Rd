\name{dell}
\Rdversion{1.1}
\alias{dell}
\alias{dellplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate coordinates of a data ellipse
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{dell} to calculates the coordinates of a 2D  data ellipse
(concentration ellipse) from (X, Y) variables.

\code{dellplus} can produce, in addition to the points of an ellipse, the
conjugate axes corresponding to a \code{chol} or other decomposition
and the surrounding parallelogram defined by these axes.

}
\usage{
dell(x, y, radius = 1, ...)

dellplus(x, y, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
Either a two-column matrix or
numeric vectors of the same length
}

  \item{radius}{
Radius of the ellipse-generating unit circle.  The default, \code{radius=1}
corresponds to a "standard" ellipse.
}
  \item{\dots}{
Other arguments passed down to \code{ell} or \code{ellplus}.
}

}

\details{
These functions simply calculate the mean vector and
covariance matrix and call \code{ell} or \code{ellplus}.
}

\value{
Returns a 2-column matrix of (X,Y) coordinates suitable for drawing with
\code{lines()}.  

For \code{dellplus}, when more than one of the options \code{ellipse}, \code{diameters},
and \code{box} is \code{TRUE}, the different parts are separated by a row of \code{NA}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Georges Monette
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cell}}, \code{\link{ell}}, \code{\link{ellplus}},
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, radius = 1, ...) 
{
    if ((is.matrix(x) && (ncol(x) > 1)) || is.data.frame(x)) 
        mat <- as.matrix(x[, 1:2])
    else if (is.list(x)) 
        mat <- cbind(x$x, x$y)
    else mat <- cbind(x, y)
    ell(apply(mat, 2, mean), var(mat), radius = radius, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{aplot}% __ONLY ONE__ keyword per line
