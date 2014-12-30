library("glmgen")

set.seed(1)
nvals=c(100, 1000, 1e4, 1e5)
kvals=c(0,1,2,3)
#nvals=c(1000)
#kvals=c(2)
nlambda=5

plot_fit = function(n,k, nlambda) {

  x = 1:n
  f = sin(x/(max(x)-min(x))*3.1*pi) 
  p = 1/(1+exp(-f))
  y = rbinom(n,1,p)
  
  fit = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=10, family="logistic",
  control=list(obj_tol=0, rho=1))
  
  par(mfrow=c(1,nlambda))
  for(i in 1:nlambda) {
    plot(x,f,type="l",col=1,
      xlab="x",ylab="f",ylim=range( c(f,fit$beta[,i]), na.rm=TRUE))
    lines(x,fit$beta[,i],col=i+1)
    
    title(paste("n=", n, ", k=", k, sprintf(",lam=%g", fit$lambda[i])))
  }
  fit  
}

pdf("plots/logi_fits.pdf", height=4, width=4*nlambda)
for(k in kvals) {
  for(n in nvals) {
    cat("k=",k,", n=", n, " ")
    time0 = proc.time()
    tryCatch({    
      plot_fit(n, k, nlambda);
    }, error = function(err) { 
      cat("Failed for n=",n,", k=",k, "\nERROR:"); print(err)       
    }, finally = {
    })
    print(proc.time()-time0)
  }
}

dev.off()

#n = 1000
#x = 1:n
#f = sin(x/(max(x)-min(x))*3.1*pi) 
#p = 1/(1+exp(-f))
#y = rbinom(n,1,p)

#nlam = 5
#time0 = proc.time()
#a0 = trendfilter(y, x, k=0, nlambda=nlam, maxiter=10, family="logistic",
#  control=list(obj_tol=0, rho=1))
#proc.time()-time0

#time1 = proc.time()
#a1 = trendfilter(y, x, k=1,  nlambda=nlam, maxiter=10, family="logistic",
#  control=list(obj_tol=0, rho=1))
#proc.time()-time0

#time2 = proc.time()
#a2 = trendfilter(y, x, k=2, nlambda=nlam, maxiter=50, family="logistic",
#  control=list(obj_tol=0, rho=1))
#proc.time()-time0


##par(ask=TRUE)
#par(mfrow=c(1,3))
#for (l in 1:nlam) {
#plot(1:nrow(a0$obj),a0$obj[,l],type="l",lty=1,main=sprintf("k=0\nlam=%0.3f",a0$lambda[l]))
#plot(1:nrow(a1$obj),a1$obj[,l],type="l",lty=1,main=sprintf("k=1\nlam=%0.3f",a1$lambda[l]))
#plot(1:nrow(a2$obj),a2$obj[,l],type="l",lty=1,main=sprintf("k=2\nlam=%0.3f",a2$lambda[l]))
#}

#pdf(file="plot_logi.pdf", width=8, height=8)
#par(ask=FALSE)
#par(mfrow=c(1,1))
##plot(x,f,type="l")
##lines(x,a1$beta[,5],col=2)
##lines(x,a2$beta[,5],col=4)

#plot(x,a2$beta[,3],type="l",col=4)
#lines(x,f,col=1)

#dev.off()

#u = exp(f)
#y = rpois(n,u)

