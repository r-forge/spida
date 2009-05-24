\name{tab.data.frame}
\Rdversion{1.1}
\alias{tab.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tab.data.frame(dd, fmla, total.margins = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dd}{
%%     ~~Describe \code{dd} here~~
}
  \item{fmla}{
%%     ~~Describe \code{fmla} here~~
}
  \item{total.margins}{
%%     ~~Describe \code{total.margins} here~~
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
function (dd, fmla, total.margins = TRUE) 
{
    if (missing(fmla)) 
        return(do.call("tab", as.list(dd, total.margins = total.margins)))
    xx = model.frame(fmla, dd, na.action = na.include)
    xx = c(xx, total.margins = total.margins)
    do.call("tab", xx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
