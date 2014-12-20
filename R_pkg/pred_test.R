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

set.seed(1)

nlambda = 3
maxiter = 100
kvals = c(0, 1, 2, 3)
nvals = c(1e3, 3e3, 1e4, 1e5)
#kvals = c(2)
even_spacing = F
thinning = 1
lam_num = 2    
x_cond = 1e10

prefix = "~/Dropbox/dropbox/code/veeru/plots/pred"

if( even_spacing ) {
  pdf(sprintf("%s_even.pdf", prefix), height=4, width=4*length(nvals))
} else pdf( sprintf("%s_uneven.pdf",prefix), height=4, width=4*length(nvals))

#k=2; n = 1000;
for(k in kvals) {
  par(mfrow=c(1,length(nvals)))
  for(n in nvals) {    
    cat("n=",n, " k=", k, "\n")
    if( even_spacing ) {  
      x = seq(0,1,length.out=n)
    } else 
      x = sort(runif(n,0,1))  
    y = sin(x*3*pi/( max(x) - min(x))) + rnorm(n,sd=0.1)

    w = rep(1,n)

    np = n
    xpred = sort(runif(np, min(x), max(x)))
    
    time0 = proc.time()    
    admm = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=maxiter, lambda.min.ratio = 1e-5, xpred=xpred, control=list(obj_tol=0, thinning=thinning, x_cond=x_cond, npred=np, lam_num=lam_num))
    time0 = proc.time() - time0

    xt = admm$x; yt = admm$y; wt = admm$w; nt = length(xt)
    ypred = admm$ypred
    #plot      
    pp_lam = lam_num
    pp_fac = max(n/500, 1) # plot a subsample of 500 points
    pp = c(1, floor((1:floor(n/pp_fac)) * pp_fac), n)
    plot(x[pp], y[pp],type="l", xlab="x",ylab="y",col="yellow",
         ylim=range(y,na.rm=TRUE));

    pp_fac = max(np/500, 1)   
    pp = c(1, floor((1:floor(np/pp_fac)) * pp_fac), np)
    lines(xpred[pp], ypred[pp], col="black")
    legend("topright",legend=c("y","predicted"),
           col=c("yellow","black"),pch=21)
    title(paste("k=", k, ", n=", n, ", #lam=", pp_lam))
  }
}

dev.off()

