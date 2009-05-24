\name{vplot}
\Rdversion{1.1}
\alias{vplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vplot(mat, type = "p", new = F, pch = 16, pop = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{new}{
%%     ~~Describe \code{new} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{pop}{
%%     ~~Describe \code{pop} here~~
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
function (mat, type = "p", new = F, pch = 16, pop = 0, ...) 
{
    help <- "\nvplot    - plots the columns of a 2 x n matrix or a vector of length 2\n         - vplot adds to the current plot resizing it to include all plotted\n           objects in a 'euclidean' frame\n         - to start a new plot, use 'new = T'\n         - to remove the last element added use 'vplot(pop=1)'\n         Associated functions:\n         - vell( mean, var) generates an ellipse, default = unit circle\n         - vbox() generates a box\n         - vobj() generates a circle in a box\n         - orthog(theta) generates an orthog matrix rotating through angle theta\n         - orthog.proj generates the matrix of an orthog. projection into span (x)\n         - vmat( .... ) generates a 2 by n matrix\n         Examples:\n           vplot( new = T )\n           vplot( vell(), 'l' )\n           vplot( cbind(c(3,1),c(1,4)) \%*\% vell())\n           vplot( pop = 1)\n           vplot( cbind(c(3,1),c(1,4)) \%*\% vell(), type = 'l', col = 'red')\n"
    if (new || !exists(".vplot")) 
        assign(".vplot", list(list(x = 0, y = 0, type = "n")), 
            pos = 1)
    a <- .vplot
    if (!missing(mat)) {
        mat <- cbind(mat)
        if (type == "v") {
            zz <- rbind(0 * mat, mat, mat/0)
            mat <- matrix(c(zz), nrow = 2)
            type = "b"
        }
        d <- dim(mat)
        if (d[1] != 2 && d[2] == 2) {
            mat <- t(mat)
            warning("mat is n x 2 and has been transposed")
        }
        a <- c(a, list(list(x = mat[1, ], y = mat[2, ], type = type, 
            pch = pch, ...)))
    }
    dat <- NULL
    for (i in seq(along = a)) {
        dat <- c(dat, a[[i]]$x, a[[i]]$y)
    }
    par(pty = "s")
    plot(range(na.omit(dat)), range(na.omit(dat)), type = "n", 
        xlab = "", ylab = "")
    if (pop > 0) {
        keep <- 1:max(1, (length(a) - (pop + 1)))
        a <- a[keep]
    }
    abline(h = 0, v = 0)
    for (i in seq(along = a)) do.call("points", a[[i]])
    assign(".vplot", a, pos = 1)
    invisible(a)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
