\name{glmgen-class}
\Rdversion{1.1}
\docType{class}
\alias{glmgen-class}
\alias{print,glmgen-method}
\alias{show,glmgen-method}
\alias{summary.glmgen-class}
\alias{print,summary.glmgen-method}
\alias{show,summary.glmgen-method}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{show-methods}
\alias{show,ANY-method}
\alias{show,classGeneratorFunction-method}
\alias{show,classRepresentation-method}
\alias{show,envRefClass-method}
\alias{show,genericFunction-method}
\alias{show,genericFunctionWithTrace-method}
\alias{show,glmgen-method}
\alias{show,MethodDefinition-method}
\alias{show,MethodDefinitionWithTrace-method}
\alias{show,MethodSelectionReport-method}
\alias{show,MethodWithNext-method}
\alias{show,MethodWithNextWithTrace-method}
\alias{show,namedList-method}
\alias{show,ObjectsWithPackage-method}
\alias{show,oldClass-method}
\alias{show,refClassRepresentation-method}
\alias{show,refMethodDef-method}
\alias{show,refObjectGenerator-method}
\alias{show,signature-method}
\alias{show,sourceEnvironment-method}
\alias{show,summary.glmgen-method}
\alias{show,traceable-method}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,glmgen-method}
\alias{print,summary.glmgen-method}

\title{Class \code{"glmgen"}}
\description{
Parent class for all estimated generalized lasso fits in
the \pkg{glmgen} package.
}
\section{Objects from the Class}{
Objects will generally be constructed as a result of a call
to estimation algorithms such as \code{\link{trendfilter}}.
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{the input observed numeric values}
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
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "glmgen")}: estimates the coefficents
    of beta for a vector of (possibly new) \code{lambda}'s. New values can
    only be estimated within the range of the previously fitted \code{lambda}
    values.}
    \item{show}{\code{signature(object = "glmgen")}: prints a summary of the
    \code{glmgen} object.
    }
	}
}
\author{
Taylor B. Arnold, Ryan Tibshirani, Veerun Sadhanala
}
\seealso{
\code{\linkS4class{trendfilter}},
}
\keyword{classes}
