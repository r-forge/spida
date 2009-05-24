\name{Lform}
\Rdversion{1.1}
\alias{Lform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lform(fit, form, data = model.frame(fit))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (fit, form, data = model.frame(fit)) 
{
    help <- "\n      Creates an L matrix using formulas evaluated in 'data' for\n      each column of the L matrix\n      Example:\n\n      library(car)\n      fit <- lm( income ~ (education + I(education^2) )* type, Prestige)\n      summary(fit)\n\n      . . .\n      Coefficients:\n                               Estimate Std. Error t value Pr(>|t|)\n       (Intercept)                891.3    23889.1   0.037  0.97032\n       education                  210.0     5638.8   0.037  0.97037\n       I(education^2)              38.3      328.3   0.117  0.90740\n       typeprof                191523.2    63022.0   3.039  0.00312 **\n       typewc                   25692.3    73888.0   0.348  0.72887\n       education:typeprof      -28133.0    10236.0  -2.748  0.00725 **\n       education:typewc         -4485.4    14007.9  -0.320  0.74956\n       I(education^2):typeprof   1017.5      451.8   2.252  0.02679 *\n       I(education^2):typewc      170.9      671.6   0.255  0.79967\n      . . .\n      \n      # estimate the marginal value of occupation for each occupation in the data set\n      \n      L <- list( 'marginal value of education' =Lform( fit,\n             form = list( 0,1,2*education, 0,0, type == 'prof', type == 'wc',\n                2*education *(type=='prof'), 2* education * (type == 'wc')),\n             data = Prestige))\n      wald( fit, L )\n      chat <- coef( wald( fit, L ), se = 2)\n      xyplot( coef +coefp+coefm ~ education | type, cbind(Prestige,chat)[order(Prestige$education),],\n              type = 'l')\n      xyplot( chat~ education | type, Prestige)\n\n"
    gg <- getFix(fit)
    Lsub <- do.call(cbind, eval(substitute(form), data))
    if ((nrow(Lsub) != nrow(data))) {
        if ((nrow(Lsub) == 1)) 
            Lsub <- Lsub[rep(1, nrow(data)), ]
        else stop("nrow(Lsub) != nrow(data)")
    }
    if (is.null(colnames(Lsub))) 
        colnames(Lsub) <- rep("", ncol(Lsub))
    L <- matrix(0, nrow = nrow(Lsub), ncol = length(gg$fixed))
    rownames(L) <- rownames(data)
    colnames(L) <- names(gg$fixed)
    Lpos <- Lsub[, colnames(Lsub) == "", drop = F]
    Lnamed <- Lsub[, colnames(Lsub) != "", drop = F]
    for (ip in seq_len(ncol(Lpos))) L[, ip] <- Lpos[, ip]
    if (ncol(Lnamed) > 0) {
        if (length(unknown <- setdiff(colnames(Lnamed), colnames(L)))) {
            stop(paste("Unknown effect(s):", unknown, collapse = " "))
        }
        for (nn in colnames(Lnamed)) L[, nn] <- Lnamed[, nn]
    }
    attr(L, "data") <- data
    L
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
