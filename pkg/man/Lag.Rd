\name{Lag}
\Rdversion{1.1}
\alias{Lag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag with respect to time within id
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lag(x, id = rep(1, length(x)), time = 1:length(x), lag = 1, at = NULL, check = TRUE, idx)
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
  \item{at}{
%%     ~~Describe \code{at} here~~
}
  \item{check}{
%%     ~~Describe \code{check} here~~
}
  \item{idx}{
%%     ~~Describe \code{idx} here~~
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
function (x, id = rep(1, length(x)), time = 1:length(x), lag = 1, 
    at = NULL, check = T, idx) 
{
    if (check) {
        comb <- paste(id, time)
        if (any(duplicated(comb))) 
            stop("id not unique in each level of time")
    }
    if (!missing(idx)) 
        time = idx
    ret <- x
    names(ret) <- paste(id, time, sep = "::")
    retnn <- if (is.null(at)) 
        paste(id, time - lag, sep = "::")
    else paste(id, at, sep = "::")
    ret[retnn]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
