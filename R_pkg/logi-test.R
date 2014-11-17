library("glmgen")

n = 100
x = 1:n
f = sin(x/(max(x)-min(x))*3.1*pi)
p = 1/(1+exp(-f))
y = rbinom(n,1,p)

nlam = 5
time0 = proc.time()
a0 = trendfilter(y, x, k=0, nlambda=nlam, maxiter=10, family="logistic",
  control=list(obj_tol=0, rho=1))
proc.time()-time0

time1 = proc.time()
a1 = trendfilter(y, x, k=1,  nlambda=nlam, maxiter=10, family="logistic",
  control=list(obj_tol=0, rho=1))
proc.time()-time0

time2 = proc.time()
a2 = trendfilter(y, x, k=2, nlambda=nlam, maxiter=10, family="logistic",
  control=list(obj_tol=0, rho=1))
proc.time()-time0

return ("done")
# par(ask=TRUE)
par(mfrow=c(1,3))
for (l in 1:nlam) {
plot(1:nrow(a0$obj),a0$obj[,l],type="l",lty=1,main=sprintf("k=0\nlam=%0.3f",a0$lambda[l]))
plot(1:nrow(a1$obj),a1$obj[,l],type="l",lty=1,main=sprintf("k=1\nlam=%0.3f",a1$lambda[l]))
plot(1:nrow(a2$obj),a2$obj[,l],type="l",lty=1,main=sprintf("k=2\nlam=%0.3f",a2$lambda[l]))
}

par(ask=FALSE)
par(mfrow=c(1,1))
plot(x,f,type="l")
lines(x,a1$beta[,5],col=2)
lines(x,a2$beta[,5],col=4)
