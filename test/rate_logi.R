library("glmgen")

set.seed(1)
nvals=c(1e3, 1e4, 3e4)
kvals=c(0, 1, 2)
nlambda=4
thinning = 1

max_iter = 50;
plot_fit = function(n,k, nlambda) {

  x = 1:n
  f = sin(x/(max(x)-min(x))*3.1*pi) 
  p = 1/(1+exp(-f))
  y = rbinom(n,1,p)

  fit = trendfilter(y, x, k=k, nlambda=nlambda, max_iter=max_iter, family="logistic", control=list(obj_tol=0, rho=1))
  
  par(mfrow=c(1,nlambda))
  for(i in 1:nlambda) {
   
    iters = 1:(fit$iter[i] - 1)
    opt_val = fit$obj[fit$iter[i],i]
    diffs = (fit$obj[iters,i] - opt_val)/ opt_val
    diffs[diffs <= 0]  = 1e-17
    plot( iters, diffs, type="l",col=1, log="y", 
      xlab="iterations",ylab="criterion rel diff")
    
    title(paste("n=", n, ", k=", k, sprintf(",lam=%d", i)))
  }
  fit  
}
pdf("plots/logi_test2.pdf", height=4, width=4*nlambda)
for(n in nvals) {
  for(k in kvals) {
    cat("k=",k,", n=", n, " ")
    time0 = proc.time()
    tryCatch({    
      fit = plot_fit(n, k, nlambda);
    }, error = function(err) { 
      cat("Failed for n=",n,", k=",k, "\nERROR:"); print(err)       
    }, finally = {
    })
    print(proc.time()-time0)
  }
}

dev.off()

