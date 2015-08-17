#' Fit a trend filtering model
#'
#' Find the trend filtering solution of some degree \code{k} for
#' an arbitrary set of penalty values \code{lambda}. Can handle
#' link functions of Gaussian, binomial, and Poisson penalized
#' loss functions.
#'
#' @param x
#'   vector of observed data locations, or when \code{y} is NULL, vector of
#'   observed responses.
#' @param y
#'   vector of observed reponses. If missing or NULL, the responses are assumed
#'   to be given through \code{x}, and the locations are assumed to be 1 through
#'   the length of \code{x}.
#' @param weights
#'   optional vector of sample weights. If missing, the weights will be assumed
#'   to be constant (unity) across all samples.
#' @param k
#'   the polynomial order of the trendfilter fit; a nonnegative integer (orders
#'   larger than 3 are not recommended). For instance, constant trend filtering
#'   (i.e., the fused lasso) uses \code{k} equal to 0, linear trend filtering uses
#'   \code{k} equal to 1, quadratic trend filtering uses \code{k} equal to 2, etc.
#' @param family
#'   the family for the link function in the trend filtering estimator. Can be
#'   either "gaussian", "logistic", or "poisson".
#' @param lambda
#'   a sequence of lambda values at which to produce a fit. Can be left blank
#'   (highly recommended for general use), at which point the algorithm will
#'   determine appropriate lambda values.
#' @param nlambda
#'   if \code{lambda} is missing, this determines the number of lambda values
#'   dynamically constructed by the algorithm.
#' @param lambda.min.ratio
#'   if \code{lambda} is missing, this determines the ratio between the largest
#'   and smallest \code{lambda} values. The values are evenly spaced on a log scale,
#'   so this ratio should typically be set fairly small.
#' @param thinning
#'   logical. If true, then the data are preprocessed so that a smaller, better
#'   conditioned data set is used for fitting. When set to \code{NULL}, the
#'   default, function will auto detect whether thinning should be applied
#'   (i.e., cases in which the numerical fitting algorithm will struggle to converge).
#' @param method
#'   the method used to calculate the fit. Currently only 'admm' is supported.
#' @param verbose
#'   logical. Should the function print out intermediate results as it is running.
#' @param control
#'   an optional named list of control parameters to pass to the underlying algorithm;
#'   see Details for more information. Names not matching any valid parameters
#'   will be silently ignored.
#'
#' @details
#'   Further algorithmic parameters can be passed by using
#'   \code{\link{trendfilter.control.list}}.
#'
#' @references
#'   Tibshirani, R. J. (2014), "Adaptive piecewise polynomial estimation
#'     via trend filtering", Annals of Statistics 42 (1): 285--323.
#'
#'   Ramdas, A. and Tibshirani R. J. (2014), "Fast and flexible ADMM algorithms
#'     for trend filtering", arXiv: 1406.2082.
#'
#' @return an object of class 'trendfilter'.
#' @author Taylor Arnold, Aaditya Ramdas, Veeranjaneyulu Sadhanala, Ryan Tibshirani
#' @seealso \code{\link{trendfilter.control.list}}
#' @useDynLib glmgen thin_R tf_R
#'
#' @examples
#'  set.seed(0)
#'  n = 100
#'  x = runif(n, min=-2*pi, max=2*pi)
#'  y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
#'  out = trendfilter(x, y, k=2)
#'
#'  xx = seq(min(x),max(x),length=100)
#'  lambda = out$lambda[25]
#'  yy = predict(out,x.new=xx,lambda=lambda)
#'  plot(x,y)
#'  lines(xx,yy,col=2)
#'
#' @export
trendfilter = function(x, y, weights, k = 2L,
                       family = c("gaussian", "logistic", "poisson"),
                       method = c("admm"),
                       lambda, nlambda = 50L, lambda.min.ratio = 1e-5,
                       thinning = NULL, verbose = FALSE,
                       control = trendfilter.control.list()) {

  cl = match.call()
  family = match.arg(family)
  method = match.arg(method)
  family_cd = match(family, c("gaussian", "logistic", "poisson")) - 1L
  method_cd = match(method, c("admm")) - 1L

  if (missing(x) || is.null(x)) stop("x must be passed.")
  if (missing(y) || is.null(y)) { y = x; x = 1L:n }
  else if (length(x) != length(y)) stop("x and y must have the same length.")
  n = length(y)
  ## orig_x = x
  ## orig_y = y
  ord = order(x)
  y = y[ord]
  x = x[ord]

  if (missing(weights)) weights = rep(1L,length(y))
  if (any(weights==0)) stop("Cannot pass zero weights.")
  ## orig_w = weights
  weights = weights[ord]

  if (is.na(family_cd)) stop("family argument must be one of 'gaussian', 'logistic', or 'poisson'.")
  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (k >= 3) warning("Large k leads to generally worse conditioning; k=0,1,2 are the most stable choices.")

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
      sX = as.double(x),
      sY = as.double(y),
      sW = as.double(weights),
      sN = length(y),
      sK = as.integer(k),
      sControl = control,
      PACKAGE = "glmgen")
    x = z$x
    y = z$y
    weights = z$w
    n = z$n    
  }

  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (missing(lambda)) {
    if (nlambda < 1L || nlambda != floor(nlambda)) stop("nlambda must be a positive integer.")
    if (lambda.min.ratio <= 0 || lambda.min.ratio >= 1) stop("lamba.min.ratio must be between 0 and 1.")
    lambda = rep(0, nlambda)
    lambda_flag = FALSE
  } else {
    if (length(lambda) == 0L) stop("Must specify at least one lambda value.")
    if (min(lambda) < 0L) stop("All specified lambda values must be nonnegative.")
    nlambda = length(lambda)
    lambda_flag = TRUE
  }
  if (!is.list(control) || (is.null(names(control)) && length(control) != 0L))
    stop("control must be a named list.")
  control = lapply(control, function(v) ifelse(is.numeric(v),
                   as.double(v[[1]]), stop("Elements of control must be numeric.")))

  z = .Call("tf_R",
    sX = as.double(x),
    sY = as.double(y),
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
    PACKAGE = "glmgen")

  if (is.null(z)) stop("Unspecified error in C code.")
  colnames(z$beta) = as.character(round(z$lambda, 3))
  
  ## if (is.null(z$obj)) z$obj = NA_real_ ## Why is this here??
## # Put back the order in which x,y,weights were given. 
## # When thinning reduces the number of points, do not put back the order
##   beta = z$beta
##   if( n == length(ord) ){ 
##     iord = order(ord)
##     y = y[iord]
##     x = x[iord]
##     weights = weights[iord]
##     beta = matrix(beta[iord,], nrow=n)
##   } else {
##   # get beta by prediction
##     beta = .Call("tf_predict_R",
##       sX = as.double(x),
##       sBeta = as.double(z$beta),
##       sN = length(x),
##       sK = as.integer(k),
##       sX0 = as.double(orig_x),
##       sN0 = length(orig_x),
##       sNLambda = length(lambda),
##       sFamily = family_cd,
##       PACKAGE = "glmgen")
##     beta = matrix(beta, nrow=n)
##   }

  out = structure(list(y = y, x = x, weights = weights, k = as.integer(k),
    lambda = z$lambda, df = z$df, beta = z$beta, family = family,
    method = method, n = length(y), p = length(y),
    m = length(y) - as.integer(k) - 1L, obj = z$obj, 
    status = z$status, iter = z$iter, call = cl),
    class = c("trendfilter","glmgen"))
  out
}

#' Control list for tuning trend filtering algorithm
#'
#' Constructs the control parameters for the trend filtering
#' algorithm. Allows the user to customize as many or as
#' little as desired.
#'
#' @param rho
#'  this is a scaling factor for the augmented Lagrangian parameter in the ADMM
#'  algorithm. To solve a given trend filtering problem with locations \code{x}
#'  at a tuning parameter value \code{lambda}, the augmented Lagrangian parameter
#'  is set to be \code{rho * lambda * ((max(x)-min(x))/n)^k}.
#' @param obj_tol
#'  the tolerance used in the stopping criterion; when the relative change in
#'  objective values is less than this value, the algorithm terminates.
#' @param max_iter
#'  number of ADMM iterations used; ignored for k=0.
#' @param max_iter_newton
#'  for non-Gaussian GLM losses, the number of outer iterations used in Newton's method.
#' @param x_cond
#'  condition number to control the degree of thinning, when applicable. Lower numbers
#'  enforce more thinning.
#' @param alpha_ls
#'  tuning parameter for the line search used in the proximal Newton procedure for
#'  non-Gaussian GLM losses.
#' @param gamma_ls
#'  tuning parameter for the line search used in the proximal Newton for non-Gaussian
#'  GLM losses.
#' @param max_iter_ls
#'  tuning parameter for the number of line search iterations in the proximal Newton
#'  procedure for non-Gaussian GLM losses.
#'
#' @return a list of parameters.
#' @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
#' @seealso \code{\link{trendfilter}}
#'
#' @examples
#'  set.seed(0)
#'  n = 100
#'  x = runif(n, min=-2*pi, max=2*pi)
#'  y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
#'  out = trendfilter(x, y, k=2, control=trendfilter.control.list(rho=3))
#'
#' @export
trendfilter.control.list = function(rho=1, obj_tol=1e-4, max_iter=200L,
                                    max_iter_newton=30L, x_cond=1e11,
                                    alpha_ls=0.5, gamma_ls=0.8, max_iter_ls=20L) {

  z <- list(rho=rho, obj_tol=obj_tol, max_iter=max_iter,
            max_iter_newton=max_iter_newton, x_cond=x_cond,
            alpha_ls=alpha_ls, gamma_ls=gamma_ls,
            max_iter_ls=max_iter_ls)
  z
}

