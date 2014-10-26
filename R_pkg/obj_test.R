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
x = seq(0,n-1,length.out=n)
y = sin(x*3*pi/n) + rnorm(n,sd=0.1)
k = 3
nlambda = 5
maxiter = 100

D = getDtfSparse(n,k)

objective = function(y, x, lam) {
  0.5 * sum((y-x)^2) + lam * sum(abs(D %*% x))
}
objectives = function(y, beta, lam) {
  nlam = length(lam);
  obj = numeric(nlam);
  for(i in 1:nlam) {
    obj[i] = objective(y, beta[,i], lam[i]); 
  }
  obj
}

u = qr.solve(t(D), y, tol=1e-14)
beta_max = y - t(D) %*% u
lambda_max = max(abs(u))
sum(abs(D %*% beta_max))
objective(y,beta_max, lambda_max)

## Path
out_path = trendpath(y, x, ord=k, maxsteps=nlambda)
lambda = out_path$lambda
obj_path = objectives(y, out_path$beta, lambda);

## New ADMM
out_admm = trendfilter(y, x, k=k, lambda=out_path$lambda, maxiter=maxiter, control=list(obj_tol=0))
obj_admm = objectives(y, out_admm$beta, lambda);

## Old ADMM
out_admm_old = trendfilterSim::trendfilter(y=y, k=k, eabs=0, erel=0, lambda=out_path$lambda, maxiter=maxiter)
obj_admm_old = objectives(y, out_admm_old$beta, lambda);

print("Objectives(NEW)")
print(obj_admm)
print("Objectives(OLD)")
print(obj_admm_old)
print("Objectives(path)")
print(obj_path)


