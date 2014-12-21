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
trendfilter(y, x, k = 0L, family = c("gaussian", "logistic", "poisson"),
            lambda, nlambda = 100L, lambda.min.ratio = 1e-05,
            thinning = TRUE, method = c("admm"),
            objective = FALSE, control = list())
}
\arguments{
  \item{y}{
vector of observed data
}
  \item{x}{
optional vector of observed data locations, which should be sorted in
ascending order. If missing, will be assumed to be the integers 1 through
the length of \code{y}.
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
very large data set is reduced to a smaller, better conditioned data.
}
  \item{method}{
the method used to calculate the fit. Currently only 'admm' is supported.
}
  \item{objective}{
logical flag indicating whether objective values of the fit should be
returned.
}
  \item{control}{
an optional named list of control parameters to pass to the underlying algorithm;
see Details for more information. Names not matching any valid parameters
will be silently ignored.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class \code{\linkS4class{trendfilter}}
}
\references{
%% ~put references to the literature/web site here ~
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