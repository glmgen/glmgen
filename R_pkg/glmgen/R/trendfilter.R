trendfilter = function(y, x, weights, k = 2L,
                        family = c("gaussian", "logistic", "poisson"),
                        method = c("admm"),
                        lambda, nlambda = 50L, lambda.min.ratio = 1e-05,
                        thinning = !missing(x), objective = TRUE,
                        verbose = FALSE, control = list()) {

  cl = match.call()
  n = length(y)
  family = match.arg(family)
  method = match.arg(method)
  family_cd = match(family, c("gaussian", "logistic", "poisson")) - 1L
  method_cd = match(method, c("admm")) - 1L

  if (!missing(x)) {
    ord = order(x)
    y = y[ord]
    x = x[ord]
  }
  if (missing(x)) x = 1L:length(y)
  if (missing(weights)) weights = rep(1L,length(y))
  if (any(weights==0)) stop("Cannot pass zero weights.")
  if (is.na(family_cd)) stop("family argument must be one of 'gaussian', 'logistic', or 'poisson'.")

  cond = (1/n) * ( (max(x) - min(x)) / min(diff(x)))^(k+1)
  if( !thinning && cond > 1e12 ) {
    warning("The x values are ill-conditioned. Consider thinning. \nSee ?trendfilter for more info.")
  }
  if (any(is.infinite(cond)) && !thinning) {
    stop("Cannot pass duplicate x values (without thinning).")
  }

  # Thin the input data:
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
  if (!objective) z$obj = matrix(NA_real_)

  out = new("trendfilter", y = y, x = x, w = weights, k = as.integer(k), lambda = z$lambda,
            beta = z$beta, family = family, method = method, n = length(y),
            p = length(y), m = length(y) - as.integer(k) - 1L, obj = z$obj,
            call = cl)
  out
}

