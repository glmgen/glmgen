library("glmgen")

n = 100
x = 1:n
f = sin(x/(max(x)-min(x))*3*pi)
u = exp(f)
y = rpois(n,u)

nlam = 20
time0 = proc.time()
a0 = trendfilter(y, x, k=0, nlambda=nlam, maxiter=20, family="poisson",
  control=list(obj_tol=0, rho=1))
proc.time()-time0

time1 = proc.time()
a1 = trendfilter(y, x, k=1, nlambda=nlam, maxiter=20, family="poisson",
  control=list(obj_tol=0, rho=1))
proc.time()-time0

time2 = proc.time()
a2 = trendfilter(y, x, k=2, nlambda=nlam, maxiter=20, family="poisson",
  control=list(obj_tol=0, rho=1))
proc.time()-time0

par(ask=TRUE)
for (l in 1:nlam) {
par(mfrow=c(1,3))
plot(1:nrow(a0$obj),a0$obj[,l],type="l",lty=1,main=sprintf("k=1\nlam=%0.3f",a1$lambda[l]))
plot(1:nrow(a1$obj),a1$obj[,l],type="l",lty=1,main=sprintf("k=1\nlam=%0.3f",a1$lambda[l]))
plot(1:nrow(a2$obj),a2$obj[,l],type="l",lty=1,main=sprintf("k=1\nlam=%0.3f",a1$lambda[l]))
}

plot(x,f,type="l")
lines(x,a1$beta[,10],col=2)
lines(x,a2$beta[,10],col=4)
