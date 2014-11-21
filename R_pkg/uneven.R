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

compare_crit = function(n,k, nlambda) {
  eps = rnorm(n,sd=0.1)
  #x = seq(0,1,length.out=n) #sort(1:n + rnorm(n,sd=0.2))#seq(0,2*pi,length.out=n)  
  #cond = mean(diff(x))/min(diff(x))
  x = sort(runif(n,0,1))  
  dx = diff(x)
  delta = 0 #mean(dx)/500;
  for(i in 1:(n-1)) {
    x[i+1] = x[i] + dx[i] + delta
  }
  cond = mean(dx+delta)/min(dx+delta)
  
  y = sin(x/(max(x)-min(x))*3*pi) + eps

  npathsteps = min(200, round( n * 0.9 ))
  maxiter = max(100, min(200, round(n*0.1)))
  
  # Path
  path = trendpath(y, x, ord=k, maxsteps=npathsteps)
  cat("\nlambda_min(path) = ", min(path$lambda), "lambda_max=", max(path$lambda),"cond=", cond, "\n")
  
  lambda.min.ratio = min(1e-5, 0.9* min(path$lambda)/max(path$lambda))

  # New ADMM
  #admm = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=maxiter, lambda.min.ratio = lambda.min.ratio, control=list(obj_tol=1e-10, rho=1) )
  admm = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=maxiter, lambda.min.ratio = lambda.min.ratio )

  crit_path = crit_admm = numeric(length(admm$lambda))

  for (i in 1:length(admm$lambda)) {
    if (admm$lambda[i] >= min(path$lambda)) {
      crit_path[i] = objective(x,y,k,admm$lambda[i],coef(path,lambda=admm$lambda[i])$beta)
    } else {
      crit_path[i] = NA
    }
    crit_admm[i] = objective(x,y,k,admm$lambda[i],admm$beta[,i])
  }

  criterions = cbind(k, n, 1:nlambda, admm$lambda, rep(NA, nlambda), rep(NA, nlambda))
  colnames(criterions)[3:6] = c("step", "lambda", "path", "admm")
  for (i in 1:nlambda) {
    criterions[i,5] = crit_path[i]
    criterions[i,6] = crit_admm[i]
  }

  criterions
#  return(list(criterions=criterions, lambda=admm$lambda))
}

set.seed(405)
kvals = c(0,1,2,3)
nvals = c(20, 100, 1e3, 1e4, 2e4)
kvals = c(2)
nvals = c(1e4)
nlambda = 50

criterions = matrix(NA, ncol=6, nrow=0)
for(k in kvals) {
  for(n in nvals) {
    cat("k=",k,", n=", n, " ")
    time0 = proc.time()
    tryCatch({    
      crit_nk = compare_crit(n, k, nlambda);
    }, error = function(err) { 
      cat("Failed for n=",n,", k=",k, "\nERROR:"); print(err) 
      crit_nk = cbind(k, n, 1:nlambda, rep(NA, nlambda), rep(NA, nlambda), rep(NA, nlambda))      
    }, finally = {
    })
    print(proc.time()-time0)

    criterions = rbind(criterions, crit_nk)
  }
}
save( criterions, nvals, kvals, nlambda, file="uneven_test.RData")

# plot
load("uneven_test.RData");
pdf("uneven_test.pdf", height=4, width=4*length(kvals))

for(i in 1:length(nvals)) {
  n = nvals[i]
  # One row for each n
  par(mfrow=c(1,length(kvals)))

  for(j in 1:length(kvals)) {
    k = kvals[j]
    offset = ((j-1) * length(nvals) + (i-1)) * nlambda 
    rows = (offset+1): (offset+nlambda)
    tryCatch({
      plot(criterions[rows,4],criterions[rows,5],log="xy",type="l",
       xlab="Lambda",ylab="Criterion",
       ylim=range(criterions[rows,5:6],na.rm=TRUE))
      points(criterions[rows,4],criterions[rows,6],col="red")
      legend("topleft",legend=c("Path","New ADMM"),
          col=c("black","red"),pch=21)
       
      title(paste("n=", n, ", k=", k))
    
    }, error = function(err) { cat("Failed for n=",n,", k=",k, "\nERROR:"); print(err) },
    finally = {    
    })  
  }
}

dev.off();

