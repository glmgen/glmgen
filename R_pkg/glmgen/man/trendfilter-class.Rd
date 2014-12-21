\name{trendfilter-class}
\Rdversion{1.1}
\docType{class}
\alias{trendfilter-class}
\alias{predict,trendfilter-method}
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
    \item{\code{family}:}{character value giving the glm family
    used for estimation. Currently supports 'gaussian', 'binomial',
    and 'poisson'.}
    \item{\code{method}:}{method used to fit the generalized lasso}
    \item{\code{n}:}{an integer giving the length of \code{y}}
    \item{\code{p}:}{an integer giving the number of columns
    in the penalty matrix \code{D}}
    \item{\code{m}:}{an integer giving the number of rows
    in the penalty matrix \code{D}}
    \item{\code{obj}:}{a, possibly missing, matrix giving the
    objective values of the fit, with one column for each value
    of lambda}
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
%%  ~~put references to the literature/web site here~~
}
\author{
Taylor B. Arnold, Ryan Tibshirani, Veerun Sadhanala
}
\seealso{
\code{\linkS4class{glmgen}},
}
\examples{
x = runif(25, min=-2*pi, max=2*pi)
y = 1.5*sin(x) + sin(2*x) + rnorm(100, sd=0.2)
out = trendfilter(y, pos, k=3)
}
\keyword{classes}
