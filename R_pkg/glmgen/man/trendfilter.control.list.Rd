\name{trendfilter.control.list}
\alias{trendfilter.control.list}
\title{
Control list for tuning trend filtering algorithm
}
\description{
Constructs the control parameters for the trend filtering
algorithm. Allows the user to customize as many or as
little as desired.
}
\usage{
trendfilter.control.list(rho=1, obj_tol=1e-6, maxiter=200L,
                          max_iter_newton=50L, x_cond=1e11,
                          alpha_ls=0.5, gamma_ls=0.8, max_iter_ls=20L)
}
\arguments{
  \item{rho}{
the standard tuning parameter in the ADMM algorithm. Scaled so that \code{rho} equal
to 1 is the default.
}
  \item{obj_tol}{
the tolerance used in the stopping criterion.
}
  \item{maxiter}{
number of admm iterations used; ignored for k=0.
}
  \item{max_iter_newton}{
for non-Gaussian GLM losses, the number of outer iterations used in Newton's method.
}
  \item{x_cond}{
condition number to control the degree of thinning, when applicable. Lower numbers
give more thinning.
}
  \item{alpha_ls}{
tuning parameter for the line search used inside of the ADMM procedure for non-Gaussian
GLM losses.
}
  \item{gamma_ls}{
tuning parameter for the line search used inside of the ADMM procedure for non-Gaussian
GLM losses.
}
  \item{max_iter_ls}{
tuning parameter for the number of iterations to use in the line search used inside of
the ADMM procedure for non-Gaussian GLM losses.
}
}
\value{
an object of class \code{\linkS4class{trendfilter}}
}
\author{
Taylor B. Arnold
}
\seealso{
\code{\linkS4class{trendfilter}}
}
\examples{
  set.seed(0)
  n = 100
  x = runif(n, min=-2*pi, max=2*pi)
  y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
  out = trendfilter(y, x, k=2, control=trendfilter.control.list(rho=3))

}