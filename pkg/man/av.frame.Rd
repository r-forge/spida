\name{av.frame}
\Rdversion{1.1}
\alias{av.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
av.frame(model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (model, ...) 
{
    help <- "\n       av.frame( model, variable)\n       returns a data frame with model.frame(model) augmented\n       by y.res and x.res, the residuals for an added variable\n       plot\n       \n       The purpose of this function is to facilitate OLS av.plots\n       for mixed models.\n\n\n\n       Example:\n\n       library(nlme)\n       library(lattice)\n       hs <- read.csv( 'http://www.math.yorku.ca/~georges/Data/hs.csv')\n       \n       # Mixed model where ses and Sex are Level 1 and Sector is Level 2\n       \n       fit.mm <- lme( mathach ~ ses * Sex * Sector, hs, random = ~ 1+ses| school)\n\n       # for diagnostics fit an OLS model using only level 1 variables interacting\n       # with the id variable\n\n       fit.ols <- lm( mathach ~ (ses * Sex ) * factor(school), hs)\n       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, 'ses:Sex'),hs), sub = 'ses:Sex')\n       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^Sex'),hs), sub = 'Sex')\n       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^ses$|^ses:f'),hs), sub = 'ses')\n\n       \n       Note : y.res is the residual from fitting the response on\n              the model matrix for fit.ols omitting any column\n              whose names is matched (as a regular expression)\n              by 'effect'\n              x.res is the residual of the first column of the\n              model matrix that is matched by 'effect' on the\n              same matrix used for y.res.\n       Caution: To make sure that the correct columns were\n              matched, the list of matched columns that are omitted\n              is printed.\n\n"
    UseMethod("av.frame")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
