
gsp.v01 <- function( x, k, degree = rep(3, length(k)+1), smooth = rep( 2, length(k)),
      intercept = FALSE, debug = FALSE ) {
    #
    # Copied from R/GeneralSplines.R on July 16, 2008
    # gsp: V 0.3:  July 16, 2008
    #        This version excludes the intercept by adding a constraint
    #        to the L matrix thus forcing the spline to evaluate to 0 at
    #        when x = 0. This allows the inclusion of an intercept in the
    #        linear model and interactions with other variables can be
    #        modeled in the usual way.
    #
    #        'gsp' is written so the parametrization does not depend on
    #        the set of values for 'x'. This allows 'gsp' to be used for
    #        prediction with values of 'x' that are different from those
    #        in the original data. This behaviour is in contrast with that
    #        of 'poly' for example.
    #
    # Original: G. Monette, June 16, 2008 (georges@yorku.ca)
    # V 0.2:    July 6, 2008
    #
    #
    # Arguments:
    #
    #    x          values where spline is evaluated
    #
    #    k          vector of knots in order
    #
    #    degree     degree of polynomial in each interval bounded by knots
    #               'degree' must satisfy:
    #                     1. length( degree ) = length( k ) + 1
    #                     2. each degree must < 5  ( to be removed in a future version )
    #
    #    smooth     degree of smoothness across each knot
    #               i.e. the highest degree required to be smooth at each knot
    #                   0  for continuity
    #                   -1 for possible discontinuity (not yet implemented)
    #
    #    intercept  if TRUE the constant term is included in the linear space
    #               spanned by the spline. If FALSE (default) the linear term
    #               is excluded and the spline evaluates to 0 when x = 0. This
    #               allow the inclusion of factors and interactions with factors
    #               in the formula for linear models.  It also allows an
    #               overall Wald test for x (e.g. using glh)
    #
    #  Caution: Basis is not orthogonal and it is desirable for numerical
    #           reasons to centre and rescale x so that its mean is close to
    #           0 and its SD is not too large.
    #
    #  Note: max order is 4
    # fits a spline with at most a quartic order
    # the smoothness at the knots and the
    # the degree of polynomials within segments can
    # be specified.
    # The resulting basis matrix is not
    # data dependent and, thus, can be used for
    # prediction for unobserved values
    #
    # Example: Fitting a quadratic up to x = 2, a linear segment from 2 to 6 and
    #          a quadratic above 6 with continuous first derivatives (i.e. smooth)
    #          across the the knots
    #
    #   a) Define the fitting spline function using gsp
    #
    #           > sp <- function( x ) gsp( x, k = c(2,6), degree = c(2,1,2), smooth = c(1,1))
    #
    #   b) Fit a model: ( assuming data for 'y', 'x' and factor 'Sex'
    #      in data frame 'dd')
    #
    #           > fit <- lm( y ~ sp(x) * Sex, dd, na.action = na.exclude)
    #           > summary(fit)
    #           > glh( fit, 'sp' )  # overall test for contribution of 'x' through spline
    #
    #   c) Display fitted values over observed values of x
    #
    #           > dd$y.hat <- predict( fit )
    #           > dd <- dd[order(dd$x),]
    #           > xyplot( y.hat ~ x , dd, groups = Sex, type = "l",
    #                 auto.key = list( columns = 2, lines = T, points = F)
    #      # or more generally
    #           > pred <- expand.grid( x = 0:30, Sex = levels( dd$Sex))
    #           > pred$y.hat <- predict( fit, pred )
    #           > xyplot( y.hat ~ x , pred, groups = Sex, type = "l",
    #                 auto.key = list( columns = 2, lines = T, points = F)
    #
    #
    #
    # Future enhancements:  Inverse function to estimate polynomial coefficients in
    #                       each interval
    #

    #
    if ( !debug ) disp = function(x) invisible(0)

    kernel = function( L )  {
        QR = qr( t(L))
        qr.Q(QR, complete = T)  [, -seq_len( QR$rank)]
    }

    Xfull = function( x, k, degree) {
        if ( length( degree) != length(k)+1) stop("length( degree ) must = length( k ) + 1")
        Xmat = function( x , degree ) {
           do.call ( "cbind" , lapply( 0:degree, function( i ) if( i == 0) rep(1,length(x)) else x^i))
        }
        k = sort(k)
        g = cut( x, c(-Inf, k, Inf))
        Xraw = Xmat(x, max(degree))
        do.call( 'cbind',
        lapply( seq_along(degree) , function( iint ) (g ==levels(g)[iint]) *  Xraw[, 1:(degree[iint]+1),drop=F]))
    }

    X0 = Xfull( 0, k, degree)

    #  Prepare 'base' set to solve for basis transformation
    if (FALSE) {
      nk = length( k )
      extend = (max(k) - min(k) + 1) / nk
      ints = c( min(k) - extend, sort(k), max(k) + extend)
      xs = approx( ints, n = (max(degree)+1)*(nk+1)) $ y
      xs = c(-xs,xs)

      X.full = Xfull( xs, k , degree)

      disp( X.full )
    }

    # Generate L matrix:


    fs = list (
      f = function(x){ c(1,x,x^2,x^3,x^4)} ,
      f1 = function(x) { c( 0, 1, 2*x, 3*x^2, 4*x^3) } ,
      f2 = function(x) { c( 0, 0, 2, 6*x, 12*x^2)} ,
      f3 = function(x) { c(0,0,0,6, 24*x)},
      f4 = function(x) { c(0,0,0,0,24)}  )

    nints = length(k) + 1

    endcols = cumsum( degree + 1 )
    disp( endcols )
    startcols = c( 1, endcols[-length(endcols)] + 1)
    disp( startcols)
    ncols = degree + 1
    disp( ncols )

    Lsmooth = matrix( 0,
        ncol = ncol(X0),
        nrow = sum( smooth + 1))

    # disp(Lsmooth)

    irow = 0

    for ( i in seq_along(k) ) {
        disp(i)
        knot = k[i]
        iplus = startcols[i]:endcols[i]
        disp(iplus)
        iminus = startcols[i+1]:endcols[i+1]
        disp(iminus)
        for ( j in 0:smooth[i] ) {
            irow = irow + 1
            Lsmooth[ irow, iplus]  =( fs[[j+1]](knot)[1:ncols[i]] )
            Lsmooth[ irow, iminus] =( - fs[[j+1]](knot)[1:ncols[i+1]] )
        }
    }


    disp(Lsmooth)
    disp( Lsmooth %*% kernel( Lsmooth))

    ## add 0 constraint
    if ( intercept == FALSE) Lsmooth = rbind( Lsmooth, X0)

    H = kernel( Lsmooth )

    Xret = Xfull( x, k, degree) %*% H
    Xret
}


##
##
##  General polynomial splines: August 1, 2008
##
##
##


##
##   Functions
##   gsp
##
Overview = "

Xmat = function( x, degree, D = 0, signif = 3)
     design/estimation matrix for f[D](x) where f(x) is polynomial of degree degree.

Xf =  function(   x, knots, degree = 3, D = 0, right = TRUE , signif = 3)
     uses Xmat to form 'full' matrix with blocks determined by knots intervals

Cmat = function( knots, degree, smooth, intercept = 0, signif = 3)
     linear constraints

Emat = function(  knots, degree, smooth , intercept = FALSE, signif = 3)
     estimates - not necessarily a basis

basis  = function( X , coef = FALSE )
     selects linear independent subset of columns of X

spline.T = function( knots, degree, smooth, intercept = 0, signif = 3 )
     full transformation of Xf to spline basis and constraints

spline.E = function( knots, degree, smooth, intercept = 0, signif = 3 )
     transformation for spline basis (first r columns of spline.T)

gsp = function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1], degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3) {
     design matrix for specified spline

"

##
##
##
##


Xmat <- function( x, degree, D = 0, signif = 3) {

        # Return design matrix for d-th derivative

        if ( length(D) < length( x ) ) D = rep(D, length.out = length(x))
        if ( length(x) < length( D ) ) x = rep(x, length.out = length(D))
        xmat = matrix( x, nrow = length(x), ncol = degree + 1)
        # disp( d)
        # disp( x)
        expvec <- 0:degree
        coeffvec <- rep(1, degree+1)
        expmat <- NULL
        coeffmat <- NULL

        for (i in 0:max(D) ) {
            expmat <- rbind(expmat, expvec)
            coeffmat <- rbind(coeffmat, coeffvec)
            coeffvec <- coeffvec * expvec
            expvec <- ifelse(expvec > 0, expvec - 1, 0)
          }
        G = coeffmat[ D + 1, ] * xmat ^ expmat[ D + 1, ]

        xlab = signif( x, signif)
        rownames(G) = ifelse( D == 0, paste('f(', xlab ,')', sep = ''), paste( "D",D,"(",xlab,")", sep = ""))
        colnames(G) = paste( "X", 0:(ncol(G)-1), sep = "")
        G
}

#  Xmat( c(1:10,pi), 5,0:4,3)

Xf =  function(   x, knots, degree = 3, D = 0, right = TRUE , signif = 3) {

    # With the default, right == TRUE, if x is at a knot then it is included in
    # in the lower interval. With right = FALSE, it is included in the higher
    # interval. This is needed when building derivative constraints at the
    # knot

    xmat = Xmat ( x, degree, D , signif )
    k = sort(knots)
    g = cut( x, c(-Inf, k, Inf), right = right)
    do.call( 'cbind',
        lapply( seq_along(levels(g)) , function( i ) (g ==levels(g)[i]) *  xmat))
}

# Xf( c(1:10,pi), c(3,6))

# Xf( 3, c(3,6),3, 0, )
# Xf( 3, c(3,6),right = F)

Cmat = function( knots, degree, smooth, intercept = 0, signif = 3) {
      # generates contraint matrix
      dm = max(degree)

      # intercept

      cmat = NULL
      if( !is.null(intercept))  cmat = rbind( cmat, "Intercept" = Xf( 0, knots, dm, D=0 ))

      # continuity constraints

      for ( i in seq_along(knots) ) {
          k = knots[i]
          sm = smooth[i]
            if ( sm > -1 ) {  # sm = -1 corresponds to discontinuity
                dmat  =Xf( k , knots, dm, D = seq(0,sm) , F ) -   Xf( k , knots, dm, D = seq(0,sm) , T )
                rownames( dmat ) = paste( "C(",signif(k, signif),").",seq(0,sm), sep = '')
                cmat = rbind( cmat,  dmat)
            }
      }

      # reduced degree constraints

      degree = rep( degree, length.out = length(knots) + 1)
      # disp ( degree )
      for ( i in seq_along( degree)) {
          di = degree[i]

          if ( dm > di ) {
                dmat = diag( (length(knots) + 1) * (dm +1)) [  (i - 1)*(dm + 1) + 1 + seq( di+1,dm), , drop = F]
                rownames( dmat ) = paste( "I.", i,".",seq(di+1,dm),sep = '')
                #disp( cmat)
                #disp( dmat)
                cmat = rbind( cmat, dmat)
          }
      }
      rk = qr(cmat)$rank
      spline.rank = ncol(cmat) - rk
      attr(cmat,"ranks") = c(npar.full = ncol(cmat), C.n = nrow(cmat), C.rank = rk , spline.rank = spline.rank )
      attr(cmat,"d") = svd(cmat) $ d
      cmat

}

# Cmat( c(-.2,.4), c(2,3,2), c(2,2))


Emat = function(  knots, degree, smooth , intercept = FALSE, signif = 3) {
      # estimation matrix:
      # polynomial of min
      # note that my intention of using a Type II estimate is not good
      #    precisely because that means the parametrization would depend
      #    on the data which I would like to avoid
      #    BUT perhaps we can implement later
      if ( length( degree) < length(knots) + 1) degree = rep( degree, length.out =  length(knots) + 1)
      dmin = min(degree)
      if ( dmin < 1) stop("minimum degree must be at least 1")
      dmax = max(degree)
      smin = min(smooth)
      smax = max(smooth)
      imax = length(degree)

      zeroi = as.numeric(cut( 0, c(-Inf, sort( knots), Inf)))
      dzero = degree[zeroi]

      cmat = Xf( 0, knots, degree = dmax , D = seq( 1, dzero))

      if ( imax > zeroi ){
          for ( i in (zeroi+1): imax) {
            d.right = degree[ i ]
            d.left  = degree[ i-1 ]
            k = knots[ i - 1 ]
            sm = smooth[ i - 1 ]
            if (  d.right > sm ) {
                dmat  =Xf( k , knots, dmax, D = seq(sm+1,d.right) , F ) -
                          Xf( k , knots, dmax, D = seq(sm+1,d.right) , T )
                rownames( dmat ) = paste( "C(",signif(k, signif),").",seq(sm+1,d.right), sep = '')
                cmat = rbind( cmat,  dmat)

            }
          }
      }

      if ( zeroi > 1 ){
          for ( i in zeroi: 2) {
            d.right = degree[ i ]
            d.left  = degree[ i-1 ]
            k = knots[ i - 1 ]
            sm = smooth[ i - 1 ]
            if (  d.left > sm ) {
                dmat  =Xf( k , knots, dmax, D = seq(sm+1,d.left) , F ) -   Xf( k , knots, dmax, D = seq(sm+1,d.left) , T )
                rownames( dmat ) = paste( "C(",signif(k, signif),").",seq(sm+1,d.left), sep = '')
                cmat = rbind( cmat,  dmat)

            }
          }
      }
      cmat
}


basis  = function( X , coef = FALSE ) {
      # returns linear independent columns
      #
      #   X = fr(X) %*% attr(fr(X),"R")
      #
        q <- qr(X)
        sel <- q$pivot[ seq_len( q$rank)]
        ret <- X[ , sel ]
        attr(ret,"cols") <- sel
        if ( coef )attr(ret,"R") <- qr.coef( qr(ret) , mat)
        colnames(ret) <- colnames(X)[ sel ]
        ret
}

#    tt( c(-1,1,2), c(3,3,3,2), c(3,3,2))

#    tmat = rbind( c(1,1,1,1,1), c(1,0,0,0,1), c(0,1,1,1,0), c(1,2,0,0,1))
#    t(tmat)
#    basis( t(tmat))

spline.T = function( knots, degree, smooth, intercept = 0, signif = 3 ) {
      cmat = Cmat( knots, degree, smooth, intercept, signif  )  # constraint matrix
      emat = Emat( knots, degree, smooth, !is.null(intercept), signif  )  # estimation matrix
      disp( list(C= cmat, E=emat ) )
      nc = nrow( cmat )
      ne = nrow( emat )
      basisT = t( basis( cbind( t(cmat), t(emat) ) ))
      cols = attr(basisT,"cols")
      ncc = sum( cols <= nc )
      Tmat = solve( basisT )
      attr( Tmat, "nC") =  ncc
      attr( Tmat, "CE") = rbind( cmat, emat)
      attr( Tmat, "CEbasis" ) = t(basisT)
      Tmat
}
# ( x =   spline.T( c(-.2,.4), c(2,3,2), c(2,2)))


spline.E = function( knots, degree, smooth, intercept = 0, signif = 3 ) {
      cmat = Cmat( knots, degree, smooth, intercept, signif  )  # constraint matrix
      emat = Emat( knots, degree, smooth, !is.null(intercept), signif  )  # estimation matrix
      # disp( list(C= cmat, E=emat ) )
      nc = nrow( cmat )
      ne = nrow( emat )
      basisT = t( basis( cbind( t(cmat), t(emat) ) ))
      cols = attr(basisT,"cols")
      ncc = sum( cols <= nc )
      Tmat = solve( basisT )
      Tmat[, (ncc+1): ncol(Tmat)]
}

# ( x =   spline.E( c(-.2,.4), c(2,3,2), c(2,2)))

gsp = function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1], degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3) {
    degree = rep( degree, length = length(knots) + 1)
    smooth = rep( smooth, length = length(knots))
    spline.attr <-   list( knots = knots, degree = degree, smooth = smooth, intercept = intercept, signif = signif)
    if (is.null(x)) return ( spline.attr)
    ret = Xf ( x, knots, max(degree), signif = signif) %*%
        spline.E( knots, degree, smooth, intercept = intercept, signif = signif)
    attr(ret, "spline.attr") <- spline.attr
    ret
}


if (FALSE) {
       x = seq( -1,1,.01)
       matplot( x, gsp( NULL, c(-.5,.5), c(2,3,2), 1), type = 'l', lwd = 2)

       1
      #  gsp( (-5):5, c(-2,3), 3)
      #  gsp( NULL, c(-2,3))

      #  y = (0:10)^3
      #  y4 =  (0:10)^4
      #  x = (-5):5

      #  summary( lm( y4 ~ gsp2( x, c(-2,3))))

      zd = data.frame( x = -10:10 , sex = sample( c('m','f'), 21, replace = T), g = sample(c('a','b'), 21, replace=T))

      zd$y = with(zd, .01*(1 + (x>0))*x^2+ (sex=='m') + x*(x>0)*(sex=='f') + .1 * rnorm(nrow(zd)))

      td( lty = c(1,2,1,2) , col = rep(c('red','blue'), each = 2), lwd = 2)
      xyplot( y ~ x , zd, groups= sex:g, auto.key= list(lines=T), type = c('p','l','p','l','g'),
          panel = panel.superpose.2)


      sp2  = function( x) gsp2( x, 0, c(2,2), 1)
      sp1  = function( x) gsp2( x, 0, c(1,1), 0)
      sp2c = function(x)  gsp2( x, .1, c(2,2), 1)
      sp2( zd$x)
      sp1( zd$x)
      fit = lm( y ~ 1+sex + sp2(x) + (sex:sp1(x)-1), zd)
      fitc = lm( y ~ 1+sex + sp2c(x) + (sex:sp1(x)-1), zd)
      summary(fit)
      summary(fitc)
      zd$y.p = predict( fit)
      xyplot( y + y.p ~ x , zd, auto.key= list(lines=T),
            panel= panel.superpose.2, type = c('p','l','g'))

      basis(cbind(sp1( zd$x), sp2( zd$x)))
      #glh( fit, ":")



      tt = function( ...) {
        rk = function ( ... ) qr(...)$rank
        cm = Cmat( ...)
        em = Emat( ...)
        ret = list( Cmat = cm, Emat = em)
        attr( ret, 'ranks') = c( p = ncol( cm),C.n=nrow(cm), C.rank = rk(cm), E.n = nrow( em),E.rank = rk(em))
        ret
      }

} # end of FALSE


##
##
##  General polynomial splines: August 8, 2008
##
##
##                                                        ][



Overview = "

Xmat = function( x, degree, D = 0, signif = 3)
     design/estimation matrix for f[D](x) where f(x) is polynomial of degree degree.

Xf =  function(   x, knots, degree = 3, D = 0, right = TRUE , signif = 3)
     uses Xmat to form 'full' matrix with blocks determined by knots intervals

Cmat = function( knots, degree, smooth, intercept = 0, signif = 3)
     linear constraints

Emat = function(  knots, degree, smooth , intercept = FALSE, signif = 3)
     estimates - not necessarily a basis

basis  = function( X , coef = FALSE )
     selects linear independent subset of columns of X

spline.T = function( knots, degree, smooth, intercept = 0, signif = 3 )
     full transformation of Xf to spline basis and constraints

spline.E = function( knots, degree, smooth, intercept = 0, signif = 3 )
     transformation for spline basis (first r columns of spline.T)

gsp = function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1],
        degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3) {
     design matrix for specified spline

"

##
##
##
##


Xmat = function( x, degree, D = 0, signif = 3) {

        # Return design matrix for d-th derivative

        if ( length(D) < length( x ) ) D = rep(D, length.out = length(x))
        if ( length(x) < length( D ) ) x = rep(x, length.out = length(D))
        xmat = matrix( x, nrow = length(x), ncol = degree + 1)
        # disp( d)
        # disp( x)
        expvec <- 0:degree
        coeffvec <- rep(1, degree+1)
        expmat <- NULL
        coeffmat <- NULL

        for (i in 0:max(D) ) {
            expmat <- rbind(expmat, expvec)
            coeffmat <- rbind(coeffmat, coeffvec)
            coeffvec <- coeffvec * expvec
            expvec <- ifelse(expvec > 0, expvec - 1, 0)
          }
        G = coeffmat[ D + 1, ] * xmat ^ expmat[ D + 1, ]

        xlab = signif( x, signif)
        rownames(G) = ifelse( D == 0, paste('f(', xlab ,')', sep = ''),
              paste( "D",D,"(",xlab,")", sep = ""))
        colnames(G) = paste( "X", 0:(ncol(G)-1), sep = "")
        G
}

#  Xmat( c(1:10,pi), 5,0:4,3)

Xf =  function(   x, knots, degree = 3, D = 0, right = TRUE , signif = 3) {

    # With the default, right ,== TRUE, if x is at a knot then it is included in
    # in the lower interval. With right = FALSE, it is included in the higher
    # interval. This is needed when building derivative constraints at the
    # knot

    xmat = Xmat ( x, degree, D , signif )
    k = sort(knots)
    g = cut( x, c(-Inf, k, Inf), right = right)
    do.call( 'cbind',
        lapply( seq_along(levels(g)) , function( i ) (g ==levels(g)[i]) *  xmat))
}

# Xf( c(1:10,pi), c(3,6))

# Xf( 3, c(3,6),3, 0, )
# Xf( 3, c(3,6),right = F)

Cmat = function( knots, degree, smooth, intercept = 0, signif = 3) {
      # generates constraint matrix
      dm = max(degree)

      # intercept

      cmat = NULL
      if( !is.null(intercept))  cmat = rbind( cmat, "Intercept" =
          Xf( 0, knots, dm, D=0 ))

      # continuity constraints

      for ( i in seq_along(knots) ) {
          k = knots[i]
          sm = smooth[i]
            if ( sm > -1 ) {  # sm = -1 corresponds to discontinuity
                dmat  =Xf( k , knots, dm, D = seq(0,sm) , F ) -   Xf( k ,
                      knots, dm, D = seq(0,sm) , T )
                rownames( dmat ) = paste( "C(",signif(k, signif),").",
                      seq(0,sm), sep = '')
                cmat = rbind( cmat,  dmat)
            }
      }

      # reduced degree constraints

      degree = rep( degree, length.out = length(knots) + 1)
      # disp ( degree )
      for ( i in seq_along( degree)) {
          di = degree[i]

          if ( dm > di ) {
                dmat = diag( (length(knots) + 1) * (dm +1)) [  (i - 1)*(dm + 1) +
                    1 + seq( di+1,dm), , drop = F]
                rownames( dmat ) = paste( "I.", i,".",seq(di+1,dm),sep = '')
                #disp( cmat)
                #disp( dmat)
                cmat = rbind( cmat, dmat)
          }
      }
      rk = qr(cmat)$rank
      spline.rank = ncol(cmat) - rk
      attr(cmat,"ranks") = c(npar.full = ncol(cmat), C.n = nrow(cmat),
              C.rank = rk , spline.rank = spline.rank )
      attr(cmat,"d") = svd(cmat) $ d
      cmat

}

# Cmat( c(-.2,.4), c(2,3,2), c(2,2))


Emat = function(  knots, degree, smooth , intercept = FALSE, signif = 3) {
      # estimation matrix:
      # polynomial of min
      # note that my intention of using a Type II estimate is not good
      #    precisely because that means the parametrization would depend
      #    on the data which I would like to avoid
      #    BUT perhaps we can implement later
      if ( length( degree) < length(knots) + 1) degree = rep( degree,
            length.out =  length(knots) + 1)
      dmin = min(degree)
      if ( dmin < 1) stop("minimum degree must be at least 1")
      dmax = max(degree)
      smin = min(smooth)
      smax = max(smooth)
      imax = length(degree)

      zeroi = as.numeric(cut( 0, c(-Inf, sort( knots), Inf)))
      dzero = degree[zeroi]

      cmat = Xf( 0, knots, degree = dmax , D = seq( 1, dzero))

      if ( imax > zeroi ){
          for ( i in (zeroi+1): imax) {
            d.right = degree[ i ]
            d.left  = degree[ i-1 ]
            k = knots[ i - 1 ]
            sm = smooth[ i - 1 ]
            if (  d.right > sm ) {
                dmat  =Xf( k , knots, dmax, D = seq(sm+1,d.right) , F ) -
                          Xf( k , knots, dmax, D = seq(sm+1,d.right) , T )
                rownames( dmat ) = paste( "C(",signif(k, signif),").",
                    seq(sm+1,d.right), sep = '')
                cmat = rbind( cmat,  dmat)

            }
          }
      }

      if ( zeroi > 1 ){
          for ( i in zeroi: 2) {
            d.right = degree[ i ]
            d.left  = degree[ i-1 ]
            k = knots[ i - 1 ]
            sm = smooth[ i - 1 ]
            if (  d.left > sm ) {
                dmat = Xf( k , knots, dmax, D = seq(sm+1,d.left) , F ) -
                      Xf( k , knots, dmax, D = seq(sm+1,d.left) , T )
                rownames( dmat ) = paste( "C(",signif(k, signif),").",
                      seq(sm+1,d.left), sep = '')
                cmat = rbind( cmat,  dmat)

            }
          }
      }
      cmat
}


basis  = function( X , coef = FALSE ) {
      # returns linear independent columns
      #
      #   X = fr(X) %*% attr(fr(X),"R")
      #
        q <- qr(X)
        sel <- q$pivot[ seq_len( q$rank)]
        ret <- X[ , sel ]
        attr(ret,"cols") <- sel
        if ( coef )attr(ret,"R") <- qr.coef( qr(ret) , mat)
        colnames(ret) <- colnames(X)[ sel ]
        ret
}

#    tt( c(-1,1,2), c(3,3,3,2), c(3,3,2))

#    tmat = rbind( c(1,1,1,1,1), c(1,0,0,0,1), c(0,1,1,1,0), c(1,2,0,0,1))
#    t(tmat)
#    basis(tmat)
#    basis( t(tmat))

spline.T =
function( knots, degree, smooth, intercept = 0, signif = 3 ) {
      cmat = Cmat( knots, degree, smooth, intercept, signif  )  # constraint matrix
      emat = Emat( knots, degree, smooth, !is.null(intercept), signif  )  # estimation matrix
      #disp( list(C= cmat, E=emat ) )
      nc = nrow( cmat )
      ne = nrow( emat )
      basisT = t( basis( cbind( t(cmat), t(emat) ) ))
      cols = attr(basisT,"cols")
      ncc = sum( cols <= nc )
      Tmat = solve( basisT )
      attr( Tmat, "nC") =  ncc
      attr( Tmat, "CE") = rbind( cmat, emat)
      attr( Tmat, "CEbasis" ) = t(basisT)
      Tmat
}

# ( x =   spline.T( c(-.2,.4), c(2,3,2), c(2,2)))


spline.E = function( knots, degree, smooth, intercept = 0, signif = 3 ) {
      cmat = Cmat( knots, degree, smooth, intercept, signif  )  # constraint matrix
      emat = Emat( knots, degree, smooth, !is.null(intercept), signif  )  # estimation matrix
      # disp( list(C= cmat, E=emat ) )
      nc = nrow( cmat )
      ne = nrow( emat )
      basisT = t( basis( cbind( t(cmat), t(emat) ) ))
      cols = attr(basisT,"cols")
      ncc = sum( cols <= nc )
      Tmat = solve( basisT )
      Tmat[, (ncc+1): ncol(Tmat)]
}

# ( x =   spline.T( c(-.2,.4), c(2,3,2), c(2,2)))
# ( x =   spline.E( c(-.2,.4), c(2,3,2), c(2,2)))

# ( x =   spline.E( c(-.2,.4), c(2,3,2), c(2,2)))

gsp = function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1],
    degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3) {
help = "
    gsp generates a matrix of regressors for a spline with
    knots, degree of polynomials in each interval and
    the degree of smoothness at each knot. Typically, gsp is
    used to define a function that is then used in a model equation.

    For example:
      simd <- data.frame( age = rep(1:50, 2), y = sin(2*pi*(1:100)/5)+ rnorm(100),
             G = rep( c('male','female'), c(50,50)))

      sp <- function(x) gsp( x, knots = c(10,25,40), degree = c(1,2,2,1),
                       smooth = c(1,1,1))

      fit <- lm( y ~ sp(age)*G, simd)
      xyplot( predict(fit) ~ age , simd, groups = G,type = 'l')
      summary(fit)

    Value:    the portion of a model matrix to fit a spline with
              given knots, polynomial degrees and degree of smoothness

    Arguments:
      x       the values for which each row of the spline matrix are generated
      knots   the knots, i.e. points at which the fitted polynomial may have
              a discontinuous derivative or second derivative, etc.
      degree  the degree of the polynomial to be used in each interval. Note
              that the number of intervals = number of knots + 1.
      smooth  the degree of smoothness at each knot, i.e. the highest degree
              for which the derivative must be continuous at the knot
      intercept (optional) use if the central value for the spline is different from 0
      signif  (optional) number of digits to retain in names of coefficients

    Details:
      A function to fit a cubic spline with knots at 5 and 10 would be

         sp <- function( x ) gsp( x, c(5,10), c(3,3,3), c(2,2))

      indicating that each a cubic polynomial is used in each of the three
      intervals and that the second derivative is continuous at each knot.

      A 'natural cubic spline' with linear components in each unbounded interval
      would have the form:

         sp <- function( x ) gsp( x, c(5,10), c(1,3,1), c(2,2))

      Finally, quadratic and linear splines, respectively:

         sp.quad <- function( x ) gsp( x, c(5,10), c(2,2,2), c(1,1))
         sp.lin  <- function( x ) gsp( x, c(5,10), c(1,1,1), c(0,0))

      Where the same degree is used for all intervals and knots, it suffices
      to give it once:

         sp.quad <- function( x ) gsp( x, c(5,10), 2, 1)
         sp.lin  <- function( x ) gsp( x, c(5,10), 1, 0)

      An easy way to specify a model in which a knot is dropped is to force
      a degree of continuity equal to the degree of adjoining polynomials, e.g.
      to drop the knot at 10, use:

         sp.1 <- function( x ) gsp( x, c(5,10), c(3,3,3), c(2,3))

      This is sometimes easier than laboriously rewriting the spline function
      for each null hypothesis.

      Depending on the maximal degree of the spline, the range of x
      should not be excessive to avoid numerical problems. The spline
      matrix generated is 'raw' and values of max(abs(x))^max(degree)
      may appear in the matrix.  For example, for a cubic spline, it
      might be desirable to rescale x and/or recenter x so abs(x) < 100
      if that is not already the case. Note that the knots need to be
      correspondingly transformed.

      The naming of coefficients should allow them to be easily interpreted.
      For example:

        >  zapsmall( gsp ( 0:10, c(3, 7) , c(2,3,2), c(1,1)) )

      D1(0) D2(0) C(3).2   C(3).3 C(7).2
f(0)      0   0.0    0.0  0.00000    0.0
f(1)      1   0.5    0.0  0.00000    0.0
f(2)      2   2.0    0.0  0.00000    0.0
f(3)      3   4.5    0.0  0.00000    0.0
f(4)      4   8.0    0.5  0.16667    0.0
f(5)      5  12.5    2.0  1.33333    0.0
f(6)      6  18.0    4.5  4.50000    0.0
f(7)      7  24.5    8.0 10.66667    0.0
f(8)      8  32.0   12.5 20.66667    0.5
f(9)      9  40.5   18.0 34.66667    2.0
f(10)    10  50.0   24.5 52.66667    4.5

      The coefficient for the first regressor is the first derivative
      at x = 0; for the second regressor, the second derivative at 0;
      the third, the saltus (change) in the second derivative at x = 3,
      the fourth, the saltus in the third derivative at x = 3 and, finally,
      the saltus in the second derivative at x = 7.

      Example:
          sp <- function(x) gsp ( x, c(3, 7) , c(2,3,2), c(1,1))

          zd <- data.frame( x = seq(0,10, .5), y = seq(0,10,.5)^2 + rnorm( 21))
          fit <- lm( y ~ sp( x ), zd)
          summary(fit)
          wald( fit, list( 'third derivatives' =
                cbind( 0, sc( sp, c(1,2,3,3,3,5,7,7,7,8), D=3,
                           type = c(0,0,0,1,2,0,0,1,2,0)))))

                  numDF denDF  F.value p.value
third derivatives     2    15 4.993575 0.02177

                 Estimate Std.Error  DF   t-value p-value Lower 0.95 Upper 0.95
  D3(1)          0.000000  0.000000  15  0.796146 0.43837   0.000000   0.000000
  D3(2)          0.000000  0.000000  15  0.796146 0.43837   0.000000   0.000000
  D3(3-)         0.000000  0.000000  15  0.796146 0.43837   0.000000   0.000000
  D3(3+)        -0.127739  0.567013  15 -0.225285 0.82480  -1.336299   1.080821
  D3(3+)-D3(3-) -0.127739  0.567013  15 -0.225285 0.82480  -1.336299   1.080821
  D3(5)         -0.127739  0.567013  15 -0.225285 0.82480  -1.336299   1.080821
  D3(7-)        -0.127739  0.567013  15 -0.225285 0.82480  -1.336299   1.080821
  D3(7+)         0.000000  0.000000 Inf       NaN     NaN   0.000000   0.000000
  D3(7+)-D3(7-)  0.127739  0.567013  15  0.225285 0.82480  -1.080821   1.336299
  D3(8)          0.000000  0.000000 Inf       NaN     NaN   0.000000   0.000000

Warning messages:
1: In min(dfs[x != 0]) : no non-missing arguments to min; returning Inf
2: In min(dfs[x != 0]) : no non-missing arguments to min; returning Inf

          Note that some coefficients that are 0 by design may lead to invalid
          DRs and t-values.
"
    degree = rep( degree, length = length(knots) + 1)
    smooth = rep( smooth, length = length(knots))
    spline.attr <-   list( knots = knots, degree = degree,
      smooth = smooth, intercept = intercept, signif = signif)
    if (is.null(x)) return ( spline.attr)
    ret = Xf ( x, knots, max(degree), signif = signif) %*%
        spline.E( knots, degree, smooth, intercept = intercept, signif = signif)
    attr(ret, "spline.attr") <- spline.attr
    ret
}


Proj.1 <- function( x)  x %*% ginv(x)

Proj <- function(x) {
   # projection matrix onto span(x)
     u <- svd(x,nv=0)$u
     u %*% t(u)

}

Proj.test <- function( x, fun = Proj) {
help = "
     # sample test:
     zz <- matrix( rnorm(120), ncol = 3) %*% diag(c(10000,1,.000001))
     Proj.test( zz, fun = Proj)
     Proj.test( zz, fun = Proj.1)
"
     # limited testing suggests Proj is generally
     # roughly as good as Proj.1  with errors 1e-17
     # but Proj.1 occasionally has errors ~ 1e-11
          pp <- fun(x)
          ret <- pp%*%pp - pp
          attr(ret,"maxabs") <- max(abs(ret))
          attr(ret,"maxabs")
}




#  round(basis( Proj(tmat)),3)
#  round(basis( Proj(t(tmat))),3)

gspf.1 = function( x, knots, degree= 3, smooth = pmax(pmin( degree[-1],
        degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3) {
    # generates a 'full' spline matrix
    tmat = spline.T(knots, degree, smooth, intercept = intercept, signif = signif)
    # put E first then C  and intercept last
    nC = attr( tmat, "nC")
    nE = ncol(tmat) - nC
    tmat = tmat[, c(seq(nC+1, ncol(tmat)), 2:nC, 1)]  # first nE columns are E matrix
    # tmat = tmat[ , - ncol(tmat)]   # drop intercept columns
    # full matrix
    X = Xf ( x, knots, max(degree), signif = signif) %*% tmat
    qqr = qr(X)
    Q = qr.Q( qqr )
    R = qr.R( qqr )
    # we need to form the linear combination of

    gsp( seq(-2,10), 0, c(2,2), 0)
    ret
}

sc <-
function( sp, x , D = 0, type = 1) {
help ='
     sc generates a portion of a hypothesis matrix for the coefficients
     of a general spline constructed with "gsp"

     sp    is the spline function for which coefficients are required
     x     values at which spline is evaluated
     D     order of derivative: 0 = value of spline, 1 = first derivative, etc.
     type  at knots: 0 limit from the left, 1 from the right, 2 is saltus (i.e. jump
           from left to right)

     Warning: sc will not work correctly if the function defining the spline
              transforms the variable, e.g. sp <- function(x) gsp( x/100, knot=2 )

     Example:

     simd <- data.frame( age = rep(1:50, 2), y = sin(2*pi*(1:100)/5)+ rnorm(100),
             G = rep( c("male","female"), c(50,50)))

     sp <- function(x) gsp( x, knots = c(10,25,40), degree = c(1,2,2,1),
                       smooth = c(1,1,1))

     fit <- lm( y ~ sp(age)*G, simd)
     xyplot( predict(fit) ~ age , simd, groups = G,type = "l")
     summary(fit)  # convenient display

# output:

Call:
lm(formula = y ~ sp(age) * G, data = simd)

Residuals:
    Min      1Q  Median      3Q     Max
-2.5249 -0.7765 -0.0760  0.7882  2.6265

Coefficients:
                      Estimate Std. Error t value Pr(>|t|)
(Intercept)           0.733267   0.605086   1.212    0.229
sp(age)D1(0)         -0.084219   0.055163  -1.527    0.130
sp(age)C(10).2        0.010984   0.006910   1.590    0.115
sp(age)C(25).2       -0.023034   0.012881  -1.788    0.077 .
Gmale                -0.307665   0.855721  -0.360    0.720
sp(age)D1(0):Gmale    0.058384   0.078012   0.748    0.456
sp(age)C(10).2:Gmale -0.010556   0.009773  -1.080    0.283
sp(age)C(25).2:Gmale  0.026410   0.018216   1.450    0.150
---


Residual standard error: 1.224 on 92 degrees of freedom
Multiple R-squared: 0.0814,     Adjusted R-squared: 0.0115
F-statistic: 1.165 on 7 and 92 DF,  p-value: 0.3308
# end of output

     L0 <- list( "hat" =
        rbind( "females at age=20" = c( 1, sc(sp,20), 0, 0* sc(sp,20)),
               "males at age=20" = c( 1, sc(sp,20), 1, 1* sc(sp,20))),
        "male-female" = rbind( "at 20" = c( 0 , 0*sc(sp,20), 1, 1*sc(sp,20))))
     wald( fit, L0 )

     ...

     L1 <- list("D(yhat)/D(age)"=
        rbind( "female at age = 25" = c(0, sc(sp,25,1), 0, 0*sc(sp,25,1)),
               "male at x = 25" = c(0, sc(sp,25,1), 0, 1*sc(sp,25,1))))
     wald( fit, L1)
# output:
               numDF denDF  F.value p.value
D(yhat)/D(age)     2    92 1.057307 0.35157

                      Estimate Std.Error DF   t-value p-value Lower 0.95 Upper 0.95
female at age = 25  0.080544  0.056974 92  1.413694 0.16083  -0.032612   0.193700
male at x = 25     -0.019412  0.056974 92 -0.340712 0.73410  -0.132568   0.093744

'
        a = sp(NULL)
        D = rep(D, length.out = length(x))
        type = rep ( type, length.out = length(x))
        left = Xf( x, knots = a$knots, degree = max( a$degree), D = D, right = T)
        right = Xf( x, knots = a$knots, degree = max( a$degree), D = D, right = F)
        cleft = c(1,0,-1) [ match( type, c(0,1,2)) ]
        cright = c(0,1,1) [ match( type, c(0,1,2)) ]
        raw = left * cleft + right * cright
        nam = rownames( raw )
        nam = sub("^f","g", nam)
        nam0 = sub( "\\)","-)", nam)
        nam1 = sub( "\\)","+)", nam)
        nam2 = paste( nam1 ,"-",nam0,sep = '')
        rownames( raw ) =
                ifelse( match( x , a$knots, 0) > 0,
                        cbind( nam0,nam1,nam2) [ cbind( seq_along(type), type+1)],
                        ifelse( type != 2, nam, '0'))
        mod = raw %*% spline.E(
                a$knots, a$degree, a$smooth,
                intercept = a$intercept,
                signif = a$signif)
        mod
}


smspline <-
function (formula, data)
{
    if (is.vector(formula)) {
        x <- formula
    }
    else {
        if (missing(data)) {
            mf <- model.frame(formula)
        }
        else {
            mf <- model.frame(formula, data)
        }
        if (ncol(mf) != 1)
            stop("formula can have only one variable")
        x <- mf[, 1]
    }
    x.u <- sort(unique(x))
    Zx <- smspline.v(x.u)$Zs
    Zx[match(x, x.u), ]
}

smspline.v <-
function (time)
{
    t1 <- sort(unique(time))
    p <- length(t1)
    h <- diff(t1)
    h1 <- h[1:(p - 2)]
    h2 <- h[2:(p - 1)]
    Q <- matrix(0, nr = p, nc = p - 2)
    Q[cbind(1:(p - 2), 1:(p - 2))] <- 1/h1
    Q[cbind(1 + 1:(p - 2), 1:(p - 2))] <- -1/h1 - 1/h2
    Q[cbind(2 + 1:(p - 2), 1:(p - 2))] <- 1/h2
    Gs <- matrix(0, nr = p - 2, nc = p - 2)
    Gs[cbind(1:(p - 2), 1:(p - 2))] <- 1/3 * (h1 + h2)
    Gs[cbind(1 + 1:(p - 3), 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
    Gs[cbind(1:(p - 3), 1 + 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
    Gs
    Zus <- t(ginv(Q))  # orig: t(solve(t(Q) %*% Q, t(Q)))
    R <- chol(Gs, pivot = FALSE)
    tol <- max(1e-12, 1e-08 * mean(diag(R)))
    if (sum(abs(diag(R))) < tol)
        stop("singular G matrix")
    Zvs <- Zus %*% t(R)
    list(Xs = cbind(rep(1, p), t1), Zs = Zvs, Q = Q, Gs = Gs,
        R = R)
}

approx.Z <-
function (Z, oldtimes, newtimes)
{
    oldt.u <- sort(unique(oldtimes))
    if (length(oldt.u) != length(oldtimes) || any(oldt.u != oldtimes)) {
        Z <- Z[match(oldt.u, oldtimes), ]
        oldtimes <- oldt.u
    }
    apply(Z, 2, function(u, oldt, newt) {
        approx(oldt, u, xout = newt)$y
    }, oldt = oldtimes, newt = newtimes)
}


smsp = function( x, knots ) {
       # generates matrix for smoothing spline over x with knots at knots
       #  Usage:   e.g. from example in smspline
       # lme ( ....., random=list(all= pdIdent(~smsp( x, 0:100) - 1),
       #              plot= pdBlocked(list(~ days,pdIdent(~smsp( x, 0:100) - 1))),
       #              Tree = ~1))

       if( max(knots) < max(x) || min(x) < min(knots)) warning( "x beyond range of knots")
       if( length( knots ) < 4 ) warning( "knots should have length at least 4")
       sp = smspline( knots)
       approx.Z( sp, knots, x)
}

if (FALSE) {
    gspf.1 ( seq( -20,20), c(-1, 4), c(3,3,3), c(2,2))
    sp = function(x) gspf.1( x, c(-1,4), c(3,3,3), c( 1,2))

    zd = data.frame( x = x <- seq( -20,20), y = x + rnorm(length(x)))
    fit = lm( y ~ sp(x), zd)
    yh1 = predict(fit)

    fit2 = lm( y ~ Xf( x, c(-1,4) , 3), zd)
    yh2 = predict(fit2)

    matplot( zd$x, resid( lm ( Xf( x , c(-1,4), 3) ~ sp(x), zd)), type = 'l')     # works

    matplot( zd$x, resid( lm ( Xf( x , c(-1,4), 3) ~ sp(x) - 1, zd)), type = 'l')     # does not work -- good

    yh1 - yh2

    sp(zd$x)

    x = -5:5
    k = c(-1, 2,3)
    Xf( x, k, 3)
    xx =Xf( x, k, 3) %*% (zs <-spline.T(k, 3, c(1,1,1)))
    dim(xx)
}



if (FALSE) {

       x = seq( -1,1,.01)
       matplot( x, gsp( NULL, c(-.5,.5), c(2,3,2), 1), type = 'l', lwd = 2)

       1
      #  gsp( (-5):5, c(-2,3), 3)
      #  gsp( NULL, c(-2,3))

      #  y = (0:10)^3
      #  y4 =  (0:10)^4
      #  x = (-5):5

      #  summary( lm( y4 ~ gsp2( x, c(-2,3))))

      zd = data.frame( x = -10:10 , sex = sample( c('m','f'), 21,
          replace = T), g = sample(c('a','b'), 21, replace=T))

      zd$y = with(zd, .01*(1 + (x>0))*x^2+ (sex=='m') + x*(x>0)*(sex=='f') + .1 * rnorm(nrow(zd)))

      td( lty = c(1,2,1,2) , col = rep(c('red','blue'), each = 2), lwd = 2)
      xyplot( y ~ x , zd, groups= sex:g, auto.key= list(lines=T), type = c('p','l','p','l','g'),
          panel = panel.superpose.2)


      sp2  = function( x) gsp( x, 0, c(2,2), 1)
      sp1  = function( x) gsp( x, 0, c(1,1), 0)
      sp2c = function(x)  gsp( x, .1, c(2,2), 1)
      sp2( zd$x)
      sp1( zd$x)
      fit = lm( y ~ 1+sex + sp2(x) + (sex:sp1(x)-1), zd)
      fitc = lm( y ~ 1+sex + sp2c(x) + (sex:sp1(x)-1), zd)
      summary(fit)
      summary(fitc)
      zd$y.p = predict( fit)
      xyplot( y + y.p ~ x , zd, auto.key= list(lines=T),
            panel= panel.superpose.2, type = c('p','l','g'))

      basis(cbind(sp1( zd$x), sp2( zd$x)))
      #glh( fit, ":")



      tt = function( ...) {
        rk = function ( ... ) qr(...)$rank
        cm = Cmat( ...)
        em = Emat( ...)
        ret = list( Cmat = cm, Emat = em)
        attr( ret, 'ranks') = c( p = ncol( cm),C.n=nrow(cm), C.rank = rk(cm), E.n = nrow( em),E.rank = rk(em))
        ret
      }

} # end of FALSE

if (FALSE ) {
    data (dapc)

    head(dapc)
    library(lattice)
    xyplot( y ~ x |sex*race, dapc, groups = id, auto.key = list(columns = 6))

    Cmat( k , c(3,3,3,3),c(2,2,2))
    ginv(spline.T( 0, c(3,3), 2))
}

