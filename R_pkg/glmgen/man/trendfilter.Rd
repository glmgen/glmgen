\name{trendfilter}
\alias{trendfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
a%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trendfilter(y, x, k = 0L, family = c("gaussian", "logistic", "poisson"), lambda, nlambda = 100L, lambda.min.ratio = 1e-05, method = c("admm", "prime_dual"), maxiter = 1e+06, objective = FALSE, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
}
  \item{lambda.min.ratio}{
%%     ~~Describe \code{lambda.min.ratio} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{objective}{
%%     ~~Describe \code{objective} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x, k, family = c("gaussian", "logistic", "poisson"),
    lambda, nlambda = 100L, lambda.min.ratio = 1e-05, method = c("admm",
        "prime_dual"), maxiter = 1e+06, objective = FALSE, control = list())
{
    cl = match.call()
    n = length(y)
    nlam = as.integer(nlambda)
    family = match.arg(family)
    method = match.arg(method)
    family_cd = match(family, c("gaussian", "logistic", "poisson")) -
        1L
    method_cd = match(method, c("admm", "prime_dual")) - 1L
    if (k < 0 || k != floor(k))
        stop("k must be a nonnegative integer.")
    if (n < k + 2)
        stop("y must have length >= k+2 for kth order trend filtering.")
    if (maxiter <= 1L)
        stop("maxiter must be greater than 1")
    if (missing(lambda)) {
        if (nlam <= 0L)
            stop("nlambda must be a positive number.")
        if (lambda.min.ratio < 0 | lambda.min.ratio > 1)
            stop("lamba.min.ratio must be between 0 and 1.")
    }
    else {
        if (length(lambda) == 0L)
            stop("Must specify at least one lambda value.")
        if (min(lambda) < 0L)
            stop("All specified lambda values must be nonnegative.")
        nlambda = length(lambda)
    }
    if (!is.list(control) | (is.null(names(control)) & length(control) !=
        0L))
        stop("control must be a named list")
    control = lapply(control, function(v) ifelse(is.numeric(v),
        as.double(v[[1]]), stop("Elements of control must be numeric.")))
    if (is.null(x))
        x = 1L:length(y)
    if (missing(lambda)) {
        lambda = rep(0, nlambda)
        lambda_flag = FALSE
    }
    else lambda_flag = TRUE
    z = .Call("tf_R", sY = as.double(y), sX = as.double(x), sN = length(y),
        sK = as.integer(k), sFamily = as.integer(family_cd),
        sMethod = as.integer(method_cd), sMaxIter = as.integer(maxiter),
        sLamFlag = as.integer(lambda_flag), sObjFlag = as.integer(objective),
        sLambda = as.double(lambda), sNlambda = as.integer(nlambda),
        sLambdaMinRatio = as.double(lambda.min.ratio), sControl = control,
        package = "glmgen")
    if (is.null(z))
        stop("Unspecified error in C code.")
    if (!is.null(z$obj))
        z$obj = z$obj[1:z$numiter]
    else z$obj = NA_real_
    colnames(z$beta) = as.character(round(z$lambda, 3))
    out = new("trendfilter", y = y, x = x, k = as.integer(k),
        lambda = z$lambda, beta = z$beta, family = family, method = method,
        n = length(y), p = length(y), d = length(y) - as.integer(k) -
            1L, obj = z$obj, call = cl)
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
