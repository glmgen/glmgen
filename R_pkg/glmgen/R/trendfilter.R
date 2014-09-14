trendfilter = function(y, x, k = 0L, family = c("gaussian", "logistic", "poisson"),
                       lambda, nlambda = 100L, lambda.min.ratio = 1e-05,
                       method = c("admm", "prime_dual"),
                       maxiter = 25L, objective = FALSE, control = list()) {

  cl = match.call()
  n = length(y)
  nlam = as.integer(nlambda)
  family = match.arg(family)
  method = match.arg(method)
  family_cd = match(family, c("gaussian", "logistic", "poisson")) - 1L
  method_cd = match(method, c("admm", "prime_dual")) - 1L

  if (missing(x)) x = 1L:length(y)
  x_cond = diff(x)
  x_cond = mean(x_cond) / min(x_cond)
  if(!is.finite(x)) stop("Cannot pass duplicate x values.")
  if(x_cond > 30)
    warning("The x values are ill-conditioned.\nSee ?trendfilter for more info.")
  if(x_cond > 100)
    error("The x values are too ill-conditioned.\nSee ?trendfilter for more info.")
  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (maxiter <= 1L) stop("maxiter must be greater than 1")
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
            sN = length(y),
            sK = as.integer(k),
            sFamily = as.integer(family_cd),
            sMethod = as.integer(method_cd),
            sMaxIter = as.integer(maxiter),
            sLamFlag = as.integer(lambda_flag),
            sObjFlag = as.integer(objective),
            sLambda = as.double(lambda),
            sNlambda = as.integer(nlambda),
            sLambdaMinRatio = as.double(lambda.min.ratio),
            sControl = control,
            package = "glmgen")

  if (is.null(z)) stop("Unspecified error in C code.")
  if (is.null(z$obj)) z$obj = NA_real_
  colnames(z$beta) = as.character(round(z$lambda, 3))

  out = new("trendfilter", y = y, x = x, k = as.integer(k), lambda = z$lambda,
            beta = z$beta, family = family, method = method, n = length(y),
            p = length(y), m = length(y) - as.integer(k) - 1L, obj = z$obj,
            call = cl)
  out
}