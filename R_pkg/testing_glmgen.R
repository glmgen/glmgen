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

objective = function(x,y, k, lambda, beta) {
  Db = diff(beta);
  if( k > 1 ) {
    for( i in 1:k) {
      Db = diff( i * (1/diff(x,i)) * Db)
    }
  }
  
  0.5 * sum((y-beta)*(y-beta)) + lambda * sum(abs(Db))
}
n = 10
k = 3
x = seq(0, 2*pi,length.out=n)
y = sin(x) + rnorm(n,sd=0.1)

out_path = trendpath(y, x, ord=k, maxsteps=25)
out_admm = trendfilter(y, x, k=k, maxiter=50, lambda=out_path$lambda)

abs_diff = abs(out_path$beta - out_admm$beta)
errs = do.call(pmax, data.frame(t(abs_diff)))
err = max(abs(out_path$beta - out_admm$beta))

print("Errors")
print(errs)
cat(sprintf("max error = %g\n", err))

# plot
colors = rainbow(3)
plot(x, y, type="n")
lines(x,y,type="l", col="black")
lines(x,out_admm$beta[,1],type="l", col="red")
lines(x,out_path$beta[,1],type="l", col="blue")

j=1
lam = out_path$lambda[j]
obj_path = objective(x,y,k,lam, out_path$beta[,j] )
obj_admm = objective(x,y,k,lam, out_admm$beta[,j] )
objs = c(obj_path, obj_admm)
print("Objectives")
print(objs)

