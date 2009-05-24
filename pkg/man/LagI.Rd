\name{LagI}
\Rdversion{1.1}
\alias{LagI}
Lag with respect to time within id with interpolation for non-integer time
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LagI(x, id, time, lag = 1, delta = 0.01, check = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{lag}{
%%     ~~Describe \code{lag} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{check}{
%%     ~~Describe \code{check} here~~
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
function (x, id, time, lag = 1, delta = 0.01, check = T) 
{
    if (check) {
        comb <- paste(id, time)
        if (any(duplicated(comb))) 
            stop("id not unique in each level of idx")
    }
    ret <- x
    id <- as.character(id)
    names(x) <- id
    names(time) <- id
    for (nn in unique(id)) {
        pos <- id == nn
        xx <- x[pos]
        tt <- time[pos]
        topred <- tt - delta
        drop <- is.na(xx) | is.na(tt)
        xxc <- xx[!drop]
        ttc <- tt[!drop]
        nl <- length(xxc)
        if (nl > 0) {
            if (nl > 1) 
                xx.p <- approx(ttc, xxc, topred)$y
            else xx.p <- NA
            xx.lag <- xx - lag * (xx - xx.p)/delta
            ret[pos] <- xx.lag
        }
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
