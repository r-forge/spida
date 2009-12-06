\name{dellplus}
\Rdversion{1.1}
\alias{dellplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data ellipse with additional features
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dellplus(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, y, ...) 
{
    if ((is.matrix(x) && (ncol(x) > 1)) || is.data.frame(x)) 
        mat <- as.matrix(x[, 1:2])
    else if (is.list(x)) 
        mat <- cbind(x$x, x$y)
    else mat <- cbind(x, y)
    ellplus(apply(mat, 2, mean), var(mat), ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
