library("genlasso")
library("trendfilterSim")
source("funs.R")
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

n = 1000
x = seq(0,1,length.out=n)
x = sort(runif(n,0,1))
dx = diff(x)
delta = mean(dx)/500;
for(i in 1:(n-1)) {
  x[i+1] = x[i] + dx[i] + delta
}
cond = mean(dx+delta)/min(dx+delta)

y = sin(x*3*pi/( max(x) - min(x))) + rnorm(n,sd=0.1)
k = 2
nlambda = 2
maxiter = 2

npathsteps = min(2, round( n * 0.9 ))
#maxiter = max(100, min(250, round(n*0.1)))

# Path
path = trendpath(y, x, ord=k, maxsteps=npathsteps)
obj_path = objectives(x, y, k, path$lambda, path$beta);
## New ADMM
admm = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=maxiter, control=list(obj_tol=0, rho=1))
obj_admm = objectives(x, y, k, admm$lambda, admm$beta);

## Old ADMM
#admm_old = trendfilterSim::trendfilter(y=y, k=k, eabs=0, erel=0, nlambda=nlambda, maxiter=maxiter)
#obj_admm_old = objectives(y, admm_old$beta, lambda);

#print("Objectives(NEW)")
print(obj_admm)

#print("Objectives(OLD)")
#print(obj_admm_old)
print("Objectives(path)")
print(obj_path)


