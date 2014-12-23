\name{glmgen-package}
\alias{glmgen-package}
\alias{glmgen}
\docType{package}
\title{
Fast generalized lasso solver
}
\description{
  An efficent implementation of several numerical algorithms
  for solving the generalized lasso. Specialized implementations are provided
  to deal with special cases, such as trendfiltering and the (graph)
  fused lasso.
}
\details{
\tabular{ll}{
Package: \tab glmgen\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2014-09-03\cr
License: \tab GPL-2 | GPL-3\cr
Depends: \tab methods\cr
}

}
\author{
Taylor B. Arnold, Ryan Tibshirani

Maintainer: Taylor Arnold <taylor.arnold@aya.yale.edu>
}
\keyword{ package }
\examples{
  x = runif(25, min=-2*pi, max=2*pi)
  y = 1.5*sin(x) + sin(2*x) + rnorm(100, sd=0.2)
  out = trendfilter(y, x, k=3)
}
