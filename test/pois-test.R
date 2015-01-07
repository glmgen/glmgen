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
  u = exp(f)
  y = rpois(n,u)
  
  fit = trendfilter(y, x, k=k, nlambda=nlambda, max_iter=10, family="poisson",
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

pdf("plots/pois_fits.pdf", height=4, width=4*nlambda)
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
