\name{predict-methods}
\docType{methods}
\alias{predict,trendfilter-method}
\alias{predict.trendfilter}
\title{
  Make function predictions given a trendfilter object
}
\description{
  Once a trendfilter object has been fit, predictions
  can be made for a new set of x points and/or a new set
  of lambda values (provided that the new lambda values lie 
  within the range of the original tuning parameter values). 
}
\usage{
\S4method{predict}{trendfilter}(object, type = c("link", "response"),
                                lambda = NULL, x.new = NULL,
                                zero_tol=1e-6)
}
\arguments{
  \item{object}{an \R object inheriting from class \code{\linkS4class{trendfilter}}.}
  \item{type}{the type of prediction required. The default is on the scale
          of the linear predictors; the alternative \code{response} is on
          the scale of the response variable.}
  \item{lambda}{option vector of lambda values at which to predict. If \code{NULL},
          the set of lambdas used for fitting are used for prediction.}
  \item{x.new}{optional vector of locations to use for prediction. If \code{NULL},
          the set of locations used for fitting are used for prediction.}
  \item{zero_tol}{tolerance for the prediction algorithm; a small number slightly
          greater than zero greatly speeds up the computations with negligable 
	  degradation to the predictions.}
}
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

  # Will throw a warning
  predict(out,x.new=-10,lambda=lambda)
}
\keyword{methods}
