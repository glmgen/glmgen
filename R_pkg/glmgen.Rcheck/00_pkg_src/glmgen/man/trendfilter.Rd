\name{trendfilter}
\alias{trendfilter}
\title{
Fit a trendfiltering model
}
\description{
Find the trendfiltering solution of some degree \code{k} for
an arbitrary set of penalty values \code{lambda}. Can handle
link functions of Gaussian, binomial, and Poisson penalized
loss functions.
}
\usage{
trendfilter(y, x, weights, k = 2L, family = c("gaussian",
            "logistic", "poisson"), method = c("admm"), lambda,
            nlambda = 50L, lambda.min.ratio = 1e-05,
            thinning = !missing(x), objective = TRUE,
            verbose = FALSE, control = trendfilter.control.list())
}
\arguments{
  \item{y}{
vector of observed data
}
  \item{x}{
optional vector of observed data locations. If missing, will be assumed to
be the integers 1 through the length of \code{y}.
}
  \item{weights}{
optional vector of sample weights. If missing, the weights will be assumed
as constant across all samples.
}
  \item{k}{
the order to the trendfilter fit; an integer larger than 0 (orders larger
than 10 are not recommended). For instance, linear trendfiltering uses \code{k}
equal to one and cubic trendfiltering uses \code{k} equal to three.
}
  \item{family}{
the family to use for the link function in the penalized estimator. Can be
either "gaussian", "logistic", or "poisson".
}
  \item{lambda}{
a sequence of lambda values at which to produce a fit. Can be left blank
(recommended for general use), at which point the algorithm will determine
appropriate lamdba values.
}
  \item{nlambda}{
if \code{lambda} is missing, this determines the number of lambda values
dynamically constructed by the algorithm.
}
  \item{lambda.min.ratio}{
if \code{lambda} is missing, this determines the ratio between the largest
and smallest \code{lambda} values. The values are distributed on a log scale,
so this ratio should typically be set fairly small.
}
  \item{thinning}{
logical. If true, the data are preprocessing so that a badly conditioned or
very large data set is reduced to a smaller, better conditioned data. Defaults
to true if and only if a set of observations are passed.
}
  \item{method}{
the method used to calculate the fit. Currently only 'admm' is supported.
}
  \item{objective}{
logical flag indicating whether objective values of the fit should be
returned.
}
  \item{verbose}{
logical. Should the function print out intermediate results as it is running.
}
  \item{control}{
an optional named list of control parameters to pass to the underlying algorithm;
see Details for more information. Names not matching any valid parameters
will be silently ignored.
}
}
\details{
Further parameters can be passed by using \code{\link{trendfilter.control.list}}.
}
\value{
an object of class \code{\linkS4class{trendfilter}}
}
\references{
  Tibshirani, R. J. and Taylor, J. (2011), "The solution path of the
  generalized lasso", Annals of Statistics 39 (3) 1335--1371.

  Tibshirani, R. J. (2014), "Adaptive piecewise polynomial estimation
  via trend filtering", Annals of Statistics 42 (1): 285--323.

  Arnold, T. B. and Tibshirani, R. J. (2014), "Efficient implementations
  of the generalized lasso dual path algorithm", arXiv: 1405.3222.
}
\author{
Taylor B. Arnold, Ryan Tibshirani, Veerun Sadhanala
}
\seealso{
\code{\linkS4class{trendfilter}}
}
\examples{
  x = runif(25, min=-2*pi, max=2*pi)
  y = 1.5*sin(x) + sin(2*x) + rnorm(100, sd=0.2)
  out = trendfilter(y, x, k=3)
}