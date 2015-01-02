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

compare_crit = function(n,k) {
  eps = rnorm(n,sd=0.1)
  x = 1:n #sort(1:n + rnorm(n,sd=0.1))#seq(0,2*pi,length.out=n)  
  y = sin(x/(max(x)-min(x))*3*pi) + eps

  nlambda = 50
  npathsteps = min(2000, round( n * 0.9 ))
  maxiter = max(100, min(500, round(n*0.1)))
  
  # Path
  path = trendpath(y, x, ord=k, maxsteps=npathsteps)
  cat("\nlambda_min(path) = ", min(path$lambda), "lambda_max=", max(path$lambda))
  
  lambda.min.ratio = min(1e-5, 0.9* min(path$lambda)/max(path$lambda))

  # New ADMM
  admm = trendfilter(y, x, k=k, nlambda=nlambda, maxiter=maxiter, control=list(obj_tol=0, rho=1), lambda.min.ratio = lambda.min.ratio)

  # Old ADMM
  admm_old = trendfilterSim::trendfilter(y=y, k=k, lambda=admm$lambda, eabs=0, erel=0, maxiter=maxiter)

  crit_path = crit_admm = crit_admm_old = numeric(length(admm$lambda))
  D = getDtfSparse(n,k)
  for (i in 1:length(admm$lambda)) {
    #crit_path[i] = objective(x,y,k,admm$lambda[i],path$beta[,i])
    if (admm$lambda[i] >= min(path$lambda)) {
      crit_path[i] = objective(x,y,k,admm$lambda[i],coef(path,lambda=admm$lambda[i])$beta)
    } else {
      crit_path[i] = NA
    }
    crit_admm[i] = objective(x,y,k,admm$lambda[i],admm$beta[,i])
    crit_admm_old[i] = objective(x,y,k,admm$lambda[i],admm_old$beta[,i])
  }

  criterions = cbind(k, n, 1:nlambda, admm$lambda, rep(NA, nlambda), rep(NA, nlambda), rep(NA, nlambda))
  colnames(criterions)[3:7] = c("step", "lambda", "path", "admm", "admm_old")
  for (i in 1:nlambda) {
    criterions[i,5] = crit_path[i]
    criterions[i,6] = crit_admm[i]
    criterions[i,7] = crit_admm_old[i]
  }

  tryCatch({
  
  par(mfrow=c(1,3))

  plot(admm$lambda,crit_path,log="xy",type="l",
       xlab="Lambda",ylab="Criterion",
       ylim=range(c(crit_path,crit_admm,crit_admm_old),na.rm=TRUE))
  points(admm$lambda,crit_admm,col="red")
  points(admm$lambda,crit_admm_old,col="blue")
  legend("topleft",legend=c("Path","New ADMM","Old ADMM"),
         col=c("black","red","blue"),pch=21)

  plot(admm$lambda,crit_admm-crit_path,col="red",
       log="x",xlab="Lambda",ylab="Criterion difference",
       ylim=range(c(crit_admm-crit_path,crit_admm_old-crit_path),na.rm=TRUE))
  points(admm$lambda,crit_admm_old-crit_path,col="blue")
  legend("topleft",legend=c("New - Path","Old - Path"),
         col=c("red","blue"),pch=21)
  title(paste("n=", n, ", k=", k))

  plot(admm$lambda,crit_admm-crit_admm_old,log="x",
       xlab="Lambda",ylab="Criterion difference")
  legend("topleft",legend=c("New - Old"),col=c("black"),pch=21)
  
  }, error = function(err) { cat("Failed for n=",n,", k=",k, "\nERROR:"); print(err) },
  finally = {    
  })
  
  criterions

}

set.seed(401)
## Test parameters
#kvals = c(0,1,2,3)
#nvals = c(20, 100, 1e3, 1e4, 2e4, 4e4, 7e4, 1e5)
kvals = c(1)
nvals = c(1e3, 1e4, 2e4)

criterions = matrix(NA, ncol=7, nrow=0)
pdf("even_n1.pdf", height=4, width=12)
for(k in kvals) {
  for(n in nvals) {
    cat("k=",k,", n=", n, " ")
    time0 = proc.time()
    crit_nk = compare_crit(n, k);
    print(proc.time()-time0)

    criterions = rbind(criterions, crit_nk)
  }
}

dev.off();
save( criterions, file="criterions_n1.RData")

criterions[!is.na(match(criterions[,3], c(1,7,14))),]
criterions[!is.na(match(criterions[,2], c(40000))),]

#par(mfrow=c(1,1))
#plot(x,y,pch=".")
#lines(x,path$beta[,5],col="black",lwd=2)
#lines(x,admm$beta[,5],col="red",lwd=2)
#legend("bottomleft",legend=c("Path","ADMM"),
#       col=c("black","red"),lty=1)


