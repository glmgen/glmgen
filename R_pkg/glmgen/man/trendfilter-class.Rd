\name{trendfilter-class}
\Rdversion{1.1}
\docType{class}
\alias{trendfilter-class}
\alias{summary,trendfilter-method}

\title{Class \code{"trendfilter"}}
\description{
Class for estimated trendfilter fits in
the \pkg{glmgen} package.
}
\section{Objects from the Class}{
Objects will generally be constructed as a result of a call
to the function \code{\link{trendfilter}}.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{the input observed numeric values}
    \item{\code{x}:}{the locations of the observed values \code{y}}
    \item{\code{w}:}{sample weights}
    \item{\code{k}:}{integer, order of the trendfilter fit}
    \item{\code{lambda}:}{vector sequency of penalty lambdas,
    sorted in decreasing order}
    \item{\code{beta}:}{matrix of beta coefficents, with one column
    for each value of \code{lambda}}
    \item{\code{family}:}{character value giving the GLM family
    used for estimation. Currently supports 'gaussian', 'binomial',
    and 'poisson'.}
    \item{\code{method}:}{method used to fit the generalized lasso}
    \item{\code{n}:}{an integer giving the length of \code{y}}
    \item{\code{m}:}{an integer giving the number of rows
    in the penalty matrix \code{D}}
    \item{\code{p}:}{an integer giving the number of columns
    in the penalty matrix \code{D}}
    \item{\code{obj}:}{a matrix giving the objective values of the fit,
    with one column for each value of lambda}
    \item{\code{call}:}{the function \code{call} used for fitting
    the generalized lasso model}
  }
}
\section{Extends}{
Class \code{"\linkS4class{glmgen}"}, directly.
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "trendfilter")}: efficent
    calculation of predictions from the fit for a new set of \code{lambda}
    and/or new set of \code{x} values. Both \code{lambda} and \code{x}
    must lie within the range of the original input values.
    }
    \item{summary}{\code{signature(object = "trendfilter")}: produces
    a summary of the trendfilter fit.
    }
  }
}
\references{
  Tibshirani, R. J. (2014), "Adaptive piecewise polynomial estimation
  via trend filtering", Annals of Statistics 42 (1): 285--323.

  Ramdas, A. and Tibshirani R. J. (2014), "Fast and flexible ADMMM algorithms
  for trend filtering", arXiv: 1406.2082.
}
\author{
Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
}
\seealso{
\code{\linkS4class{glmgen}},
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
}
\keyword{classes}
