\name{capply}
\Rdversion{1.1}
\alias{capply}
\alias{cvar}
\alias{dvar}
\alias{dvar.factor}
\alias{dvar.default}
\alias{cvar.factor}
\alias{cvar.default}
\concept{multilevel model}
\concept{mixed model}
\concept{hierarchical data}
\concept{hierarchical model}
\concept{contextual variable}
\concept{contextual effect}
\concept{compositional effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a function within each cluster of multilevel data
}
\description{
\code{capply} that can be used with vectors as with tapply
or with data.frames to produce 'contextual variables' computed
within subgroups of the data and expanded to a constant over
elements of each subgroup.
}
\usage{
capply( x, by, FUN, ...)

capply( df, form, FUN, ...)

cvar( x, by, ... )

dvar( x, by, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector or data frame that provides the first argument of \code{FUN}
%%     ~~Describe \code{x} here~~
}
  \item{by} {
            If \code{x} is a vector: a 'factor' of the same lenth as \code{x} whose
            levels identify clusters.
            If \code{x} is a data frame, a one-sided formula that identifies the variable(s)
            within \code{x} to be used to identify clusters.
  }
  \item{FUN} { a function to be applied to \code{x} within each cluster. \code{FUN} can return a
    single value, or a vector whose length is equal to the number of elements in each cluster.
  }
  \item{\dots}{  additional variables to be supplied to \code{FUN}
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Function that can be used with vectors as with tapply
or with data.frames to produce 'contextual variables' computed
within subgroups of the data and expanded to a constant over
elements of each subgroup.

\code{capply( x , by, FUN , ...)} where \code{x} is a vector
  is equivalent to
  
\code{unsplit ( lapply ( split ( x , by ), FUN, ...), by )}

  which has the same effect as

\code{tapply( x, by, FUN, ...) [ tapply( x, by) ]} if FUN returns a vector of length 1.
If \code{FUN} returns a vector, it is recycled to the length of the input value.

When the first argument is a data frame:

capply ( dd, by, FUN, ...) uses unsplit - lapply - split
  to apply \code{FUN} to each sub data frame. In this case, by can
  be a formula that is evaluated in 'dd'.
  
  This syntax makes it easy to compute formulas involving more than one
  variable in 'dd'. An example:
  
  \code{capply( dd, ~gg, function(x) with( x, mean(Var1) / mean(Var2) ) )}
  
  where 'Var1' and 'Var2' are numeric variables and 'gg' a grouping factor in
  data frame 'dd'.

\code{cvar} is intended to create contextual variable in model formulas. If 'x' is numerical
it is equivalent to \code{capply(x,id,mean)}. If 'x' is a factor, 'cvar' generates the
equivalent of a model matrix for the factor with indicators replaced by the proportion
within each cluster.

\code{dvar} is equivalent to \code{x - cvar(x,by)} and creates what is commonly known as
a version of 'x' that is 'centered within groups' (CWG).  It creates the correct matrix
for a factor so that the between group interpretation of the effect of \code{cvar(x,by)}
is that of the 'between group' or 'compositional' effect of the factor.
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

     data( hs )
     head( hs )

     # FUN returns a single value
     hs$ses.mean <- capply( hs$ses, hs$school, mean, na.rm = T)
     hs$ses.hetero <- capply ( hs$ses, hs$school, sd , na.rm = T)
     hs.summ <- up( hs, ~school )
     head( summ )   # variables invariant within school
     
     # FUN returns a vector
     # with 'x' a data frame
     # Note how the 'with' function provides an easy way to write use a
     #   formula as the '...' variable.
     
     hs$minority.prop <- capply( hs, ~ school, with, mean( Minority == "Yes"))
     
     # equivalently:
     
     hs$minority.prop <- capply( hs$Minority, hs$school, mean)
     
     # on very large data frames with many columns that are not used, the 'data frame'
     # version of 'capply' can be very slow in comparison with 'vector' version.
     
     # In contrast with 'tapply' 'FUN' can return a vector, e.g. ranks within groups
     
     hs$mathach.rank <- capply( hs, ~ school, with , rank(mathach))
     
     # cvar and dvar in multilevel models
     
     library( nlme )
     data ( hs )
     fit <- lme( mathach ~ Minority * Sector, hs, random = ~ 1 | school)
     summary ( fit )
     
     fit.contextual <- lme( mathach ~ (Minority + cvar(Minority, school)) * Sector,
                       hs, random = ~ 1| school)
     summary(fit.contextual) # contextual effect of cvar(Minority)
     
     fit.compositional <- lme( mathach ~ (dvar(Minority,school) + cvar(Minority, school)) * Sector,
                       hs, random = ~ 1| school)
     summary(fit.compositional) # compositional effect of cvar(Minority)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
