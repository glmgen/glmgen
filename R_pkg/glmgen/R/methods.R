#' Get coefficients from a glmgen object
#'
#' @method coef glmgen
#'
#' @param object
#'   output of summary.iolm
#' @param lambda
#'   optional vector of lambda values to calculate coefficients
#'   at. If missing, will use break points in the fit.
#' @param ...
#'   optional, currently unused, arguments
#'
#' @export
coef.glmgen = function (object, lambda = NULL, ...) {
  # If no lambda given, just return beta
  if (is.null(lambda))
    return(object$beta)

  # If all lambda are equal to some computed lambda, just
  #   return propely transformed version of object$beta
  if (all(!is.na(index <- match(lambda, object$lambda))))
    return(object$beta[,index,drop=FALSE])

  if (min(lambda) < 0) stop("All specified lambda values must be nonnegative.")
  if (min(lambda) < min(object$lambda) | max(lambda) > max(object$lambda))
    stop("Cannot predict lambda outside the range used when fitting.")

  # If here, need to interpolate lambdas
  o = order(lambda, decreasing = TRUE)
  o2 = order(object$lambda, decreasing=TRUE)
  lambda = lambda[o]
  knots = object$lambda[o2]
  k = length(lambda)
  mat = matrix(rep(knots, each = k), nrow = k)
  b = lambda >= mat
  blo = max.col(b, ties.method = "first")
  bhi = pmax(blo - 1, 1)
  i = bhi == blo
  p = numeric(k)
  p[i] = 0
  p[!i] = ((lambda - knots[blo])/(knots[bhi] - knots[blo]))[!i]

  betas = object$beta[,o2, drop=FALSE]
  beta = t((1 - p) * t(betas[, blo, drop = FALSE]) +
      p * t(betas[, bhi, drop = FALSE]))
  colnames(beta) = as.character(round(lambda, 3))

  return(beta[, order(o), drop=FALSE])
}

#' Get predictions from a trendfilter object
#'
#' @method predict trendfilter
#'
#' @param object
#'   output of summary.iolm
#' @param lambda
#'   optional vector of lambda values to calculate coefficients
#'   at. If missing, will use break points in the fit.
#' @param ...
#'   optional, currently unused, arguments
#'
#' @export
predict.trendfilter = function (object, lambda = NULL, ...) {
  # If no lambda given, just return beta
  if (is.null(lambda))
    return(object$beta)

  # If all lambda are equal to some computed lambda, just
  #   return propely transformed version of object$beta
  if (all(!is.na(index <- match(lambda, object$lambda))))
    return(object$beta[,index,drop=FALSE])

  if (min(lambda) < 0) stop("All specified lambda values must be nonnegative.")
  if (min(lambda) < min(object$lambda) | max(lambda) > max(object$lambda))
    stop("Cannot predict lambda outside the range used when fitting.")

  # If here, need to interpolate lambdas
  o = order(lambda, decreasing = TRUE)
  o2 = order(object$lambda, decreasing=TRUE)
  lambda = lambda[o]
  knots = object$lambda[o2]
  k = length(lambda)
  mat = matrix(rep(knots, each = k), nrow = k)
  b = lambda >= mat
  blo = max.col(b, ties.method = "first")
  bhi = pmax(blo - 1, 1)
  i = bhi == blo
  p = numeric(k)
  p[i] = 0
  p[!i] = ((lambda - knots[blo])/(knots[bhi] - knots[blo]))[!i]

  betas = object$beta[,o2, drop=FALSE]
  beta = t((1 - p) * t(betas[, blo, drop = FALSE]) +
      p * t(betas[, bhi, drop = FALSE]))
  colnames(beta) = as.character(round(lambda, 3))

  return(beta[, order(o), drop=FALSE])
}

#' Print the output of a glmgen object
#'
#' @method print glmgen
#'
#' @param x
#'   object to print
#' @param ...
#'   optional, currently unused, arguments
#'
#' @export
print.glmgen = function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  cat("\nOutput:\n")
  cat(paste0(class(x), " model with ", length(x$lambda),
      " values of lambda.", "\n\n"))
}

#' Print the output of a glmgen summary object
#'
#' @method print summary.glmgen
#'
#' @param x
#'   object to print
#' @param ...
#'   optional, currently unused, arguments
#'
#' @export
print.summary.glmgen = function(x, ...) {
  class(x) = "matrix"
  print(x, digits = 4, print.gap = 3)
}

#' Summarize a generic glmgen object
#'
#' @method summary glmgen
#'
#' @param object
#'   object to summarize
#' @param ...
#'   optional, currently unused, arguments
#'
#' @export
summary.glmgen = function(object, ...) {
  class(object) = "matrix"
  print(object, digits = 4, print.gap = 3)
}

#' Summarize a trendfilter glmgen object
#'
#' @method summary trendfilter
#'
#' @param object
#'   object to summarize
#' @param ...
#'   optional, currently unused, arguments
#'
#' @export
summary.trendfilter = function(object, ...) {
  df = apply(object$beta != 0, 2, sum)
  rss = colSums((object$y - predict(object, type="response"))^2)
  mat = cbind(df, object$lambda, rss)
  rownames(mat) = rep("", nrow(mat))
  colnames(mat) = c("df", "lambda", "rss")
  class(mat) = "summary.glmgen"
  return(mat)
}

