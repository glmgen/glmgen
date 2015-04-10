\name{glmgen-package}
\alias{glmgen-package}
\alias{glmgen}
\docType{package}
\title{
  Fast generalized lasso algorithms
}
\description{
  Efficient algorithms for generalized lasso problems. Specialized
  implementations are provided to deal with particular problem classes,  
  such as trend filtering.
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
Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani

Maintainer: Taylor Arnold <taylor.arnold@acm.org>
}
\keyword{ package }
\examples{
  set.seed(0)
  n = 100
  x = runif(n, min=-2*pi, max=2*pi)
  y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
  out = trendfilter(y, x, k=2)

  xx = seq(min(x),max(x),length=500)
  lambda = out@lambda[25]
  yy = predict(out,x.new=xx,lambda=lambda)
  plot(x,y)
  lines(xx,yy,col=2)
}
