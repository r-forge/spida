\name{ell}
\Rdversion{1.1}
\alias{ell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ell(center, shape, radius = 1, n = 100)
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
function (center, shape, radius = 1, n = 100) 
{
    fac <- function(x) {
        xx <- svd(x)
        t(xx$v) * sqrt(pmax(xx$d, 0))
    }
    angles = (0:n) * 2 * pi/n
    if (length(radius) > 1) {
        ret <- lapply(radius, function(r) rbind(r * cbind(cos(angles), 
            sin(angles)), NA))
        circle <- do.call(rbind, ret)
    }
    else circle = radius * cbind(cos(angles), sin(angles))
    ret <- t(c(center) + t(circle \%*\% fac(shape)))
    attr(ret, "parms") <- list(center = rbind(center), shape = shape, 
        radius = radius)
    class(ret) <- "ell"
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
