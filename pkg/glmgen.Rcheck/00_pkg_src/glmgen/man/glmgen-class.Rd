\name{glmgen-class}
\Rdversion{1.1}
\docType{class}
\alias{glmgen-class}
\alias{coef,glmgen-method}
\alias{print,glmgen-method}
\alias{show,glmgen-method}

\title{Class \code{"glmgen"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("glmgen", ...)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{y}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{beta}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{family}:}{Object of class \code{"character"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{p}:}{Object of class \code{"integer"} ~~ }
    \item{\code{d}:}{Object of class \code{"integer"} ~~ }
    \item{\code{obj}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "glmgen")}: ... }
    \item{print}{\code{signature(x = "glmgen")}: ... }
    \item{show}{\code{signature(object = "glmgen")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("glmgen")
}
\keyword{classes}
