\name{brace}
\Rdversion{1.1}
\alias{brace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a brace for a graph
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
brace(x1 = 0, y1 = 0, x2 = 0, y2 = 1, right = T, rad = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
}
  \item{right}{
%%     ~~Describe \code{right} here~~
}
  \item{rad}{
%%     ~~Describe \code{rad} here~~
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
function (x1 = 0, y1 = 0, x2 = 0, y2 = 1, right = T, rad = 0.2) 
{
    uin <- par("uin")
    ux <- uin[1]
    uy <- uin[2]
    dx <- x2 - x1
    dy <- y2 - y1
    alpha <- atan(ux * dx, uy * dy)
    scale <- sqrt((ux * dx)^2 + (uy * dy)^2)
    if (scale > 5 * rad) 
        rad <- rad/scale
    qcirc <- cbind(cos((0:10) * pi/20), sin((0:10) * pi/20))
    qcircr <- cbind(cos((10:0) * pi/20), sin((10:0) * pi/20))
    rot <- function(theta) t(cbind(c(cos(theta), sin(theta)), 
        c(-sin(theta), cos(theta))))
    seg1 <- t(t(rad * qcirc \%*\% rot(-pi/2)) + c(0, rad))
    seg4 <- t(t(rad * qcirc) + c(0, 1 - rad))
    seg3 <- t(t((rad * qcircr) \%*\% rot(pi)) + c(2 * rad, 0.5 + 
        rad))
    seg2 <- t(t((rad * qcircr) \%*\% rot(pi/2)) + c(2 * rad, 0.5 - 
        rad))
    bra <- rbind(seg1, seg2, seg3, seg4)
    if (!right) 
        bra <- bra \%*\% diag(c(-1, 1))
    bra <- scale * bra \%*\% rot(-alpha)
    bra <- bra \%*\% diag(c(1/ux, 1/uy))
    bra <- t(t(bra) + c(x1, y1))
    bra
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
