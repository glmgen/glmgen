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


nlambda = 3
maxiter = 100
nvals = c(1e4, 1e5, 3e5, 1e6)
kvals = c(1,2,3)
x_conds = c(1e7, 1e8, 1e9, 1e11)
nvals = c(1e4)
kvals = c(2)
x_conds = c(1e11)
#presmooth = 1 # hacky way of setting minavg thinning
thinning = 1
even_spacing = T

set.seed(1)

results = matrix(NA, ncol=8, nrow=0)
# k n C nt time step lam obj

if( even_spacing ) {
  pdf("~/Dropbox/dropbox/code/veeru/plots/thin_test.pdf", height=4, width=4*length(nvals))
} else pdf("~/Dropbox/dropbox/code/veeru/plots/thin_test.pdf", height=4, width=4*length(nvals))

#k=1; n = 2000; x_cond = 1e17
for(k in kvals) {
  par(mfrow=c(1,length(nvals)))  
  for(n in nvals) {    
    for(x_cond in x_conds) {
     cat("n=",n, " k=", k, "x_cond=", x_cond, "\n")
      if( even_spacing ) {  
        x = seq(0,1,length.out=n)
      } else 
        x = sort(runif(n,0,1))  
      y = sin(x*3*pi/( max(x) - min(x))) + rnorm(n,sd=0.1)
      
      w = rep(1,n)

      time0 = proc.time()    
      admm = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=maxiter, lambda.min.ratio = 1e-2, control=list(obj_tol=0, thinning=thinning, x_cond=x_cond))
      time0 = proc.time() - time0
      
      xt = admm$x; yt = admm$y; wt = admm$w; nt = length(xt)
      obj_admm = objectives(xt, yt, wt, k, admm$lambda, admm$beta);    
      results_nk = cbind(k, n, x_cond, nt, time0[3], 1:nlambda, admm$lambda, obj_admm)
      colnames(results_nk)[4:8] = c("nt", "time", "step", "lambda", "criterion")
      
      results = rbind(results, results_nk)
      
      #plot      
      pp_lam = 2
      pp_fac = max(n/500, 1) # plotting a subsample of 500 points
      pp = c(1, floor((1:floor(n/pp_fac)) * pp_fac), n)
      plot(x[pp], y[pp],type="l", xlab="x",ylab="y",col="yellow",
       ylim=range(y,na.rm=TRUE));
      pp_fac = max(nt/500, 1)
      pp = c(1, floor((1:floor(nt/pp_fac)) * pp_fac), nt)
      lines(xt[pp], admm$beta[pp,pp_lam], col="black")
      legend("topright",legend=c("y","beta"),
          col=c("yellow","black"),pch=21)
      title(paste("k=", k, ", n=", n, ", C=", x_cond, "nt=", nt))
         
      cat("min(dx)=", min(diff(xt)), "\n")
      print("Objectives"); print(obj_admm)

    }
  }
}

dev.off()

if( even_spacing ) {
  save( results, nvals, kvals, x_conds, nlambda, file="~/Dropbox/dropbox/code/veeru/plots/thin_test.RData")
} else
  save( results, nvals, kvals, x_conds, nlambda, file="~/Dropbox/dropbox/code/veeru/plots/thin_test.RData")

# 1. For evenly spaced x in [0,1], what is the n for which ADMM starts to break?
# k=2: 5e5, more divergence for higher lambda
# k=3: 1e5
# Rule of thumb looks like n should satisfy n^{k+1} < 1/DBL_EPSILON for a successful run

# 2. Unevenly spaced points with controlled conditioning
#   C denotes Condition num mean(dx)/min(dx)
#   C = 1 gives answers similar to the even case
#   C = 30 avoids divergence until n=5e4 x=unif(n,0,1) 
#   C = 50 avoids until n=2e4
#   C>50 etc. does not avoid divergence because the divergence is due to the very large value of lambda * D^{k+1} and C does not help reduce it by much.
#   C = inf (without conditioning) divergence happens as early as n=3000.
#   No difference if x is uniformly drawn from (0,n) instead of (0,1)
#  

# 3. Averaging


