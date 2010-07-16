\name{ellplus}
\Rdversion{1.1}
\alias{ellplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ellipse with additional features
}
\description{
Calculates coordinates of an
ellipse with additional features such as center, conjugate axes, subtending parallelogram.
}
\usage{
ellplus(center = rep(0, 2), shape = diag(2), radius = 1, n = 100, angles = (0:n) * 2 * pi/n, fac = chol, ellipse = all, diameters = all, box = all, all = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{center}{
X,Y location of the center of the ellipse
}
  \item{shape}{
A 2x2 matrix, typically a covariance matrix of data (for a data ellipse), or a
covariance matrix of estimated parameters in a model (for a confidence ellipse).
}
  \item{radius}{
Radius of the ellipse-generating unit circle
}
  \item{n}{
Number of points on the unit circle used to calculate the ellipse
}
  \item{angles}{
Angles around the unit circle used to calculate the ellipse
}
  \item{fac}{
A function defining the conjugate axes used to transform the unit circle
into an ellipse.  
}
  \item{ellipse}{
Logical to indicate if the points on the ellipse should be returned
}
  \item{diameters}{
Logical to indicate if the points defining the ends of the
conjugate axes of the ellipse should be returned
}
  \item{box}{
Logical to indicate if the points on the conjugate-axes bounding box should be returned
}
  \item{all}{
Logical to request \code{ellipse}, \code{diameters} and \code{box}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a 2-column matrix of (X,Y) coordinates suitable for drawing with
\code{lines()}.  When more than one of the options \code{ellipse}, \code{diameters},
and \code{box} is \code{TRUE}, the different parts are separated by a row of \code{NA}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Georges Monette
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cell}}, \code{\link{dell}}, \code{\link{dellplus}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (center = rep(0, 2), shape = diag(2), radius = 1, n = 100, 
    angles = (0:n) * 2 * pi/n, fac = chol, ellipse = all, diameters = all, 
    box = all, all = F) 
{
    help <- "\n        ellplus can produce, in addition to the points of an ellipse, the\n        conjugate axes corresponding to a chol or other decomposition\n        and the surrounding parallelogram.\n        "
    rbindna <- function(x, ...) {
        if (nargs() == 0) 
            return(NULL)
        if (nargs() == 1) 
            return(x)
        rbind(x, NA, rbindna(...))
    }
    if (missing(ellipse) && missing(diameters) && missing(box)) 
        all <- TRUE
    circle <- function(angle) cbind(cos(angle), sin(angle))
    Tr <- fac(shape)
    ret <- list(t(c(center) + t(radius * circle(angles) \%*\% Tr)), 
        t(c(center) + t(radius * circle(c(0, pi)) \%*\% Tr)), t(c(center) + 
            t(radius * circle(c(pi/2, 3 * pi/2)) \%*\% Tr)), t(c(center) + 
            t(radius * rbind(c(1, 1), c(-1, 1), c(-1, -1), c(1, 
                -1), c(1, 1)) \%*\% Tr)))
    do.call("rbindna", ret[c(ellipse, diameters, diameters, box)])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
