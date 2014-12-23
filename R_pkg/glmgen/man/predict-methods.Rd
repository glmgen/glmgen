\name{predict-methods}
\docType{methods}
\alias{predict,trendfilter-method}
\title{ Make predictions given a trendfilter object }
\description{
  Once a trendfilter object has been fit, predictions
  can be made for a new set of x points and a new set
  of lambdas, provided each are within the range of
  the original input values.
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
          the set of lambda used in fitting are used for prediction.}
  \item{x.new}{optional vector of locations to use for prediction. If \code{NULL},
          the set of locations used for fitting are used for prediction.}
  \item{zero_tol}{tolerance for the prediction algorithm; a small number slightly
          greater than zero greatly speds up the fitting with negligable degradation
          to the fit.}
}
\examples{
  x = 1:25
  y = 1.5*sin(x) + sin(2*x) + rnorm(100, sd=0.2)
  out = trendfilter(y, x, k=3)

  pred = predict(out)
}
\keyword{methods}