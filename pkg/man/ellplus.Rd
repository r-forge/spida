\name{ellplus}
\Rdversion{1.1}
\alias{ellplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ellipse with additional features
}
\description{
Ellipse with additional features such as center, axes, subtending parallelogram.
}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ellplus(center = rep(0, 2), shape = diag(2), radius = 1, n = 100, angles = (0:n) * 2 * pi/n, fac = chol, ellipse = all, diameters = all, box = all, all = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{center}{
%%     ~~Describe \code{center} here~~
}
  \item{shape}{
%%     ~~Describe \code{shape} here~~
}
  \item{radius}{
%%     ~~Describe \code{radius} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{angles}{
%%     ~~Describe \code{angles} here~~
}
  \item{fac}{
%%     ~~Describe \code{fac} here~~
}
  \item{ellipse}{
%%     ~~Describe \code{ellipse} here~~
}
  \item{diameters}{
%%     ~~Describe \code{diameters} here~~
}
  \item{box}{
%%     ~~Describe \code{box} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
