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
n = 6
eps = rnorm(n,sd=0.1)
x = seq(0,1,length.out=n)
y = sin(x*3*pi/( max(x) - min(x) )) 
k = 3
max_iter = 100

# Path
out_path = trendpath(y, x, ord=k, maxsteps=2)
out_admm = glmgen::trendfilter(y, x, k=k, nlambda=2, max_iter=max_iter, control=list(obj_tol=0))
out_admm_old = trendfilterSim::trendfilter(y=y, k=k, nlambda=2, eabs=0, erel=0, max_iter=max_iter)

cat(sprintf("path lambda_max = %g\n",out_path$lambda[1]))
cat(sprintf("NEW lambda_max = %g\n",out_admm$lambda[1]))
#cat(sprintf("OLD lambda_max = %g\n",out_admm_old$lambda[1]))

D = getDtfPosSparse(n,k,x)
print(D)
#beta_max = y - t(D) %*% qr.solve(DDt, D %*% y)
lambda_max = max( abs(qr.solve(t(D), y)))
cat(sprintf("lambda_max = %g\n",lambda_max))


