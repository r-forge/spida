\name{abind}
\Rdversion{1.1}
\alias{abind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
abind(arr1, arr2, d, facename = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr1}{
%%     ~~Describe \code{arr1} here~~
}
  \item{arr2}{
%%     ~~Describe \code{arr2} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{facename}{
%%     ~~Describe \code{facename} here~~
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
function (arr1, arr2, d, facename = "") 
{
    help <- "\nabind                coursefun.R     for PSYC 6140/MATH 6630 05/06\n\nBinds two conformable arrays along a dimension\n\nDescription:\n\n     'abind' binds two conformable arrays, arr1 and arr2 along\n     the dimension d.\n\nUsage:\n\n     abind( arr1, arr2 , d ,  facename )\n\nArguments:\n\n     arr1, arr2: arrays, matrices or vectors\n\n       d: dimension along which arr1 and arr2 are joined\n\n     ...: Name of extended index if arr2 has one fewer dimensions than arr1\n\nDetails:\n\n      dim( arr1 ) and dim( arr2 ) must be equal except in the dth dimension.\n      If the length of dim( arr2 ) is 1 less than that of dim( arr1 ), then\n      'arr2' is treated as if it had dropped the dth dimension with size 1.\n\nValue:\n\n     The returned value, ret, is an array with dimension dim(arr1)\n     except for the dth dimension where dim(ret)[d] == dim(arr1)[d] + dim(arr2)[d].\n     If length(dim(arr2)) == length(dim(arr1)) - 1, then arr2 is treated as if\n     it dropped a dimension of size 1 in the dth dimension. 'facename' is used\n     as the name of the dimnames list for this dimension.\n\nNotes:\n\n     abind is used by atotal\n\nReferences:\n\nContributed by:  G. Monette  2005-10-10\n\nModifications:\n\n"
    d1 <- dim(arr1)
    n1 <- length(d1)
    dn1 <- dimnames(arr1)
    d2 <- dim(arr2)
    n2 <- length(d2)
    dn2 <- dimnames(arr2)
    arenull <- is.null(dn1) & is.null(dn2)
    if (is.null(dn1)) {
        dn1 <- lapply(as.list(d1), function(x) seq(1, x))
        dimnames(arr1) <- dn1
    }
    if (n1 != n2) {
        d2 <- d1
        d2[d] <- 1
        dn2 <- dn1
        dn2[[d]] <- facename
        dimnames(arr2) <- NULL
        dim(arr2) <- d2
        dimnames(arr2) <- dn2
        n2 <- n1
    }
    if (is.null(dn2)) {
        dn2 <- lapply(as.list(d2), function(x) seq(1, x))
        dimnames(arr2) <- dn2
    }
    perm <- 1:n1
    perm[c(d, n1)] <- c(n1, d)
    arr.p1 <- aperm(arr1, perm)
    arr.p2 <- aperm(arr2, perm)
    dret <- d1[perm]
    dret[n1] <- dret[n1] + d2[d]
    dnret <- dn1
    dnret[[d]] <- c(dnret[[d]], dn2[[d]])
    ret <- c(arr.p1, arr.p2)
    dim(ret) <- dret
    ret <- aperm(ret, perm)
    dimnames(ret) <- dnret
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
