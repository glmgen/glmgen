library(glmgen)
source("funs.R")

addtf = function(y, x, ks, lambda, max_iter = 10L, verbose = FALSE) {

  n = nrow(x)
  d = ncol(x)

  beta = matrix(0,n,d)
  obj = rep(0,max_iter)

  # Backfitting
  yc = y - mean(y)

  for (t in 1:max_iter) {
    for (j in 1:d) {
      cat( "t=", t, "\tj=", j, "\n")
      fit = trendfilter( y=yc-rowSums(beta)+beta[,j], x = x[,j], k=ks[j], lambda = lambda)
      # print( dim(fit@beta) )
      # if thinning reduces # of points, then the following line fails
      beta[,j] = fit@beta
    }

    # Compute objective
    obj_t = 0.5*sum((yc - rowSums(beta))^2);
    for (j in 1:d) {
      obj_t = obj_t + lambda * sum(abs(calcDx(x=x[,j], beta=beta[,j],k=ks[j])))
    }
    obj[t] = obj_t;
  }
  if( verbose ) {
    print(obj)
  }

  out <- list(beta = beta, obj = obj)
  out
}

