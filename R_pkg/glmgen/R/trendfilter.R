#source("thinning.R")
trendfilter = function(y, x, weights, k = 2L, family = c("gaussian", "logistic", "poisson"),
                       lambda, nlambda = 50L, lambda.min.ratio = 1e-05,
                       method = c("admm", "prime_dual"),
                       maxiter = 100L, control = list()) {

  cl = match.call()
  n = length(y)
  nlam = as.integer(nlambda)
  family = match.arg(family)
  method = match.arg(method)
  family_cd = match(family, c("gaussian", "logistic", "poisson")) - 1L
  method_cd = match(method, c("admm", "prime_dual")) - 1L

  if (missing(x)) x = 1L:length(y)
  if (missing(weights)) weights = rep(1L,length(y))
  if (any(weights==0)) stop("Cannot pass zero weights.")

  cond = (1/n) * ( (max(x) - min(x)) / min(diff(x)))^(k+1)
  thinning = FALSE
  x_cond = 1e10
  if( "thinning" %in% names(control) ) thinning = control$thinning
  if( "x_cond" %in% names(control)) x_cond = control$x_cond
  
  if( !thinning && cond > x_cond ) {
    warning("The x values are ill-conditioned. Consider presmoothing. \nSee ?trendfilter for more info.")
  }
  if( cond > x_cond && thinning ) {
    thinned = thin(x,y,weights,k,x_cond)
    x = thinned$x; y = thinned$y; weights = thinned$w
  }

  if(any(!is.finite(cond))) stop("Cannot pass duplicate x values.\nUse observation weights instead.")
  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (maxiter < 1L) stop("maxiter must be a positive integer")
  if (missing(lambda)) {
    if (nlam <= 0L) stop("nlambda must be a positive number.")
    if (lambda.min.ratio < 0 | lambda.min.ratio > 1)
      stop("lamba.min.ratio must be between 0 and 1.")
    lambda = rep(0, nlambda)
    lambda_flag = FALSE
  } else {
    if (length(lambda) == 0L) stop("Must specify at least one lambda value.")
    if (min(lambda) < 0L) stop("All specified lambda values must be nonnegative.")
    nlambda = length(lambda)
    lambda_flag = TRUE
  }
  if (!is.list(control) | (is.null(names(control)) & length(control) != 0L))
    stop("control must be a named list")
  control = lapply(control, function(v) ifelse(is.numeric(v),
                   as.double(v[[1]]), stop("Elements of control must be numeric.")))

  z = .Call("tf_R",
            sY = as.double(y),
            sX = as.double(x),
            sW = as.double(weights),
            sN = length(y),
            sK = as.integer(k),
            sFamily = as.integer(family_cd),
            sMethod = as.integer(method_cd),
            sMaxIter = as.integer(maxiter),
            sLamFlag = as.integer(lambda_flag),
            sLambda = as.double(lambda),
            sNlambda = as.integer(nlambda),
            sLambdaMinRatio = as.double(lambda.min.ratio),
            sControl = control,
            package = "glmgen")

  # if (is.null(z)) stop("Unspecified error in C code.")
  # if (is.null(z$obj)) z$obj = NA_real_
  # colnames(z$beta) = as.character(round(z$lambda, 3))

  ## out = new("trendfilter", y = y, x = x, w = weights, k = as.integer(k), lambda = z$lambda,
  ##   beta = z$beta, family = family, method = method, n = length(y),
  ##   p = length(y), m = length(y) - as.integer(k) - 1L, obj = z$obj,
  ##   call = cl)
  ## out

  z[["x"]] = x
  z[["y"]] = y
  z[["w"]] = weights
  z
}

# x are assumed to be sorted
thin = function(x,y,w,k,x_cond) {
  cat("Thinning interval wise\n")
  n = length(x)
  r = x[n] - x[1]
  delta = r * (n*x_cond)^(-1/(k+1))
  if(min(diff(x)) >= delta) return(list(x=x,y=y,w=w))

  # Divide the range of x into m intervals of same size
  m = min( floor(r/delta), 5*n )
  delta = r/m
    
  if( m <= 1 ) stop("Not thinning as it merges all points into one")
    
  intvl = rep(0,n)
  for(j in 1:n) {
    intvl_xj = floor( (x[j] -x[1]) / delta ) + 1 # All intvls are of type [) except the last one
    intvl[j] = max(1, min(intvl_xj, m))
  }

  # number of intervals with at least one point
  mm = length( unique(intvl) )
  xt = rep(0,mm)
  yt = rep(0,mm)
  wt = rep(0,mm)  

  lo = 1; hi = 1
  i = 1 # range 1:mm
  cur_intvl = 1 # range 1:m
  for(j in 1:n) {  
    
    if( intvl[j] > cur_intvl ) {  # crossed the current interval

      hi = j-1      
      #if( i ==  ) cat("lo=", lo, "hi=", hi, "\n")
      wt[i] = sum( w[lo:hi] )
      #xt[i] = sum( w[lo:hi] * x[lo:hi] ) / wt[i]
      xt[i] = x[1] + (cur_intvl-0.5) * delta
      yt[i] = sum( w[lo:hi] * y[lo:hi] ) / wt[i]
      i = i+1 
      cur_intvl = intvl[j]
      lo = j
    }
    if( i >= mm ) { # last interval
      hi = n
      wt[mm] = sum( w[lo:hi] ) 
      xt[mm] = sum( w[lo:hi] * x[lo:hi] ) / wt[mm]
      yt[mm] = sum( w[lo:hi] * y[lo:hi] ) / wt[mm]      
      break;
    }
  }
  cat("min(diff(x) after thinning = ", min(diff(xt)), "\n") 
  return(list(x=xt,y=yt,w=wt))
}

x_cond_num = function(x,k) {
  return(length(x) * ( (x[n] - x[1]) / min(diff(x)) )^(k+1))
}
