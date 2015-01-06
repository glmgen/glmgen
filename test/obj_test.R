library("genlasso")
library("trendfilterSim")
library("glmgen")
source("funs.R")

trendpath = function(...) {
  x = genlasso::trendfilter(...)
  class(x) = c("trendpath", class(x)[-1])
  x
}
plot.trendpath = function(x, ...) {
  class(x) = c("trendfilter", class(x)[-1])
  genlasso::plot.trendfilter(x, ...)
}
#set.seed(100)
n = 10000
x = seq(0,1,length.out=n)
x = sort(runif(n,0,1))
y = sin(x*3*pi/( max(x) - min(x))) + rnorm(n,sd=0.1)
k = 2
nlambda = 5
max_iter = 100

npathsteps = min(nlambda, round( n * 0.9 ))
#max_iter = max(100, min(250, round(n*0.1)))

# Path
path = trendpath(y, x, ord=k, maxsteps=npathsteps)
obj_path = objectives(x, y, k, path$lambda, path$beta);
## New ADMM
admm = trendfilter(y, x, k=k, lambda=path$lambda, max_iter=max_iter, control=list(obj_tol=1e-10))
obj_admm = objectives(x, y, k, admm$lambda, admm$beta);

print("Objectives(NEW)")
print(obj_admm)

print("Objectives(path)")
print(obj_path)


