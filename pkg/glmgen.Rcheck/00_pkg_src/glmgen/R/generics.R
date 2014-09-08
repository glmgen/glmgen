# Generic functions for pkg glmgen
setGeneric("coef")
setGeneric("predict")
setGeneric("print")
setGeneric("show")
setGeneric("summary")

# # Methods for S4 classes defined in glmgen
setMethod("coef", signature(object = "glmgen"),
          function(object, lambda = NULL) {
            # If no lambda given, just return beta
            if (is.null(lambda))
              return(object@beta)

            # If all lambda are equal to some computed lambda, just
            #   return propely transformed version of object@beta
            if (all(!is.na(index <- match(lambda, object@lambda))))
              return(object@beta[,index,drop=FALSE])

            if (min(lambda) < 0) stop("All specified lambda values must be nonnegative.")
            if (min(lambda) < min(object@lambda) | max(lambda) > max(object@lambda))
              stop("Cannot predict lambda outside the range used when fitting.")

            # If here, need to interpolate lambdas
            o = order(lambda, decreasing = TRUE)
            o2 = order(object@lambda, decreasing=TRUE)
            lambda = lambda[o]
            knots = object@lambda[o2]
            k = length(lambda)
            mat = matrix(rep(knots, each = k), nrow = k)
            b = lambda >= mat
            blo = max.col(b, ties.method = "first")
            bhi = pmax(blo - 1, 1)
            i = bhi == blo
            p = numeric(k)
            p[i] = 0
            p[!i] = ((lambda - knots[blo])/(knots[bhi] - knots[blo]))[!i]

            betas = object@beta[,o2, drop=FALSE]
            beta = t((1 - p) * t(betas[, blo, drop = FALSE]) +
                p * t(betas[, bhi, drop = FALSE]))
            colnames(beta) = as.character(round(lambda, 3))

            return(beta[, order(o), drop=FALSE])
          })

setMethod("predict", signature(object = "trendfilter"),
          function(object, type = c("link", "response"), lambda = NULL, x.new = NULL) {
            type = match.arg(type)
            co = coef(object, lambda)
            if (is.null(lambda)) lambda = object@lambda
            if (!is.null(x.new)) {
              if (min(x.new) < min(object@x) | max(x.new) > max(object@x))
                stop("Cannot predict outside of the original x range.")

              o = order(x.new, decreasing = TRUE)
              o2 = order(object@x, decreasing=TRUE)
              x.new = x.new[o]
              knots = object@x[o2]
              k = length(x.new)
              mat = matrix(rep(knots, each = k), nrow = k)
              b = x.new >= mat
              blo = max.col(b, ties.method = "first")
              bhi = pmax(blo - 1, 1)
              i = bhi == blo
              p = numeric(k)
              p[i] = 0
              p[!i] = ((x.new - knots[blo])/(knots[bhi] - knots[blo]))[!i]

              xs = co[o2,,drop=FALSE]
              x = t((1 - p) * t(xs[blo,, drop = FALSE]) +
                  p * t(xs[bhi,, drop = FALSE]))
              colnames(x) = as.character(round(lambda, 3))
              rownames(x) = as.character(round(x.new, 3))

              x = x[order(o),,drop=FALSE]
            } else x = co

            if (object@family == "gaussian" | type == "link")
              return(x)
            else if (object@family == "logistic")
              return(1 / (1 + exp(-x)))
            else
              return(exp(x))
          })

setMethod("print", signature(x = "glmgen"),
          function(x) {
            cat("\nCall:\n")
            dput(x@call)
            cat("\nOutput:\n")
            cat(paste0(class(x), " model with ", length(x@lambda),
                " values of lambda.", "\n\n"))
          })

setMethod("print", signature(x = "summary.glmgen"),
          function(x) {
            class(x) = "matrix"
            print(x, digits = 4, print.gap = 3)
          })

setMethod("show", signature(object = "glmgen"),
          function(object) {
            cat("\nCall:\n")
            dput(object@call)
            cat("\nOutput:\n")
            cat(paste0(class(object), " model with ",
                length(object@lambda), " values of lambda.", "\n"))
            cat("\nFamily:\n")
            dput(object@family)
            cat("\nMethod:\n")
            dput(object@method)
            cat("\n")
          })

setMethod("show", signature(object = "summary.glmgen"),
          function(object) {
            class(object) = "matrix"
            print(object, digits = 4, print.gap = 3)
          })

setMethod("summary", signature(object = "trendfilter"),
          function(object) {
            df = apply(object@beta != 0, 2, sum)
            rss = colSums((object@y - predict(object, type="response"))^2)
            mat = cbind(df, object@lambda, rss)
            rownames(mat) = rep("", nrow(mat))
            colnames(mat) = c("df", "lambda", "rss")
            new("summary.glmgen", summary = mat)
            return(mat)
          })

