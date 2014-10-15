library("genlasso")
library("glmgen")

trendpath = function(...) {
  x = genlasso::trendfilter(...)
  class(x) = c("trendpath", class(x)[-1])
  x
}

plot.trendpath = function(x, ...) {
  class(x) = c("trendfilter", class(x)[-1])
  genlasso::plot.trendfilter(x, ...)
}

n = 10
x = seq(0,2*pi,length.out=n)
y = sin(x) + rnorm(n,sd=0.1)

out_path = trendpath(y, x, ord=1, maxsteps=25)
out_admm = trendfilter(y, x, k=1, lambda=out_path$lambda)

max(abs(out_path$beta[,1] - out_admm$beta[,1]))