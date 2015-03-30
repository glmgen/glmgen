trendfilter = function(y, x, weights, k = 2L,
                        family = c("gaussian", "logistic", "poisson"),
                        method = c("admm"),
                        lambda, nlambda = 50L, lambda.min.ratio = 1e-05,
                        thinning = NULL, verbose = FALSE,
                        control = trendfilter.control.list()) {

  cl = match.call()
  family = match.arg(family)
  method = match.arg(method)
  family_cd = match(family, c("gaussian", "logistic", "poisson")) - 1L
  method_cd = match(method, c("admm")) - 1L

  n = length(y)
  if (!missing(x) && length(x)!=n) stop("x and y must have the same length.")
  if (missing(x)) x = 1L:n
  ord = order(x)
  y = y[ord]
  x = x[ord]

  if (missing(weights)) weights = rep(1L,length(y))
  if (any(weights==0)) stop("Cannot pass zero weights.")
  weights = weights[ord]

  if (is.na(family_cd)) stop("family argument must be one of 'gaussian', 'logistic', or 'poisson'.")
  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (k > 3) warning("Large k leads to generally worse conditioning; k=0,1,2 are the most stable choices.")

  cond = (1/n) * ((max(x) - min(x)) / min(diff(x)))^(k+1)
  if (!is.null(thinning) && !thinning && any(is.infinite(cond))) {
    stop("Cannot pass duplicate x values; use observation weights, or turn on thinning.")
  }
  if( !is.null(thinning) && !thinning && cond > control$x_cond ) {
    warning("The x values are ill-conditioned. Consider thinning. \nSee ?trendfilter for more info.")
  }

  # Thin the input data:
  if (is.null(thinning)) {
    if (cond > control$x_cond) {
      thinning = TRUE
    } else {
      thinning = FALSE
    }
  }

  if (thinning) {
    z = .Call("thin_R",
          sY = as.double(y),
          sX = as.double(x),
          sW = as.double(weights),
          sN = length(y),
          sK = as.integer(k),
          sControl = control,
          package = "glmgen")
    y = z$y
    x = z$x
    weights = z$w
    n = z$n    
  }

  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (missing(lambda)) {
    if (nlambda < 1L || nlambda != floor(nlambda)) stop("nlambda must be a positive integer.")
    if (lambda.min.ratio <= 0 | lambda.min.ratio >= 1) stop("lamba.min.ratio must be between 0 and 1.")
    lambda = rep(0, nlambda)
    lambda_flag = FALSE
  } else {
    if (length(lambda) == 0L) stop("Must specify at least one lambda value.")
    if (min(lambda) < 0L) stop("All specified lambda values must be nonnegative.")
    nlambda = length(lambda)
    lambda_flag = TRUE
  }
  if (!is.list(control) | (is.null(names(control)) & length(control) != 0L))
    stop("control must be a named list.")
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
            sLamFlag = as.integer(lambda_flag),
            sLambda = as.double(lambda),
            sNlambda = as.integer(nlambda),
            sLambdaMinRatio = as.double(lambda.min.ratio),
            sVerbose = as.integer(verbose),
            sControl = control,
            package = "glmgen")

  if (is.null(z)) stop("Unspecified error in C code.")
  if (is.null(z$obj)) z$obj = NA_real_
  colnames(z$beta) = as.character(round(z$lambda, 3))

# Put back the order in which x,y,weights were given. 
# When thinning reduces the number of points, do not put back the order
  beta = z$beta
  if( n == length(ord) ){ 
    iord = order(ord)
    y = y[iord]
    x = x[iord]
    weights = weights[iord]
    beta = matrix(beta[iord,])
  }

  
  out = new("trendfilter", y = y, x = x, w = weights, k = as.integer(k),
            lambda = z$lambda, beta = beta, family = family,
            method = method, n = length(y), p = length(y),
            m = length(y) - as.integer(k) - 1L, obj = z$obj,
            status = z$status, iter = z$iter, call = cl)

  out
}

trendfilter.control.list = function(rho=1, obj_tol=1e-6, max_iter=200L,
                          max_iter_newton=50L, x_cond=1e11,
                          alpha_ls=0.5, gamma_ls=0.8, max_iter_ls=20L) {

  z <- list(rho=rho, obj_tol=obj_tol, max_iter=max_iter,
            max_iter_newton=max_iter_newton, x_cond=x_cond,
            alpha_ls=alpha_ls, gamma_ls=gamma_ls,
            max_iter_ls=max_iter_ls)
  z
}

