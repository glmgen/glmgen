library(glmgen)
source("funs.R")

compute_obj = function(y, x, ks, lambda, beta) {

  d = ncol(x)
  yc = y - mean(y)
  obj = 0.5*sum((yc - rowSums(beta))^2);
  for (j in 1:d) {
    obj = obj + lambda * sum(abs(calcDx(x=x[,j], beta=beta[,j],k=ks[j])))
  }
  obj
}

addtf = function(y, x, ks, lambda, max_iter = 10L, verbose = FALSE, obj_tol = 1e-8) {

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
      beta[,j] = fit@beta
    }

    obj[t] = compute_obj(y=y,x=x,ks=ks, lambda=lambda, beta=beta);

    if (t > 1) {
      if( abs(obj[t] - obj[t-1]) <= abs(obj[t]) * obj_tol ) {
        break
      }
    }
  }
  if( verbose ) {
    print(obj)
  }

  out <- list(beta = beta, obj = obj, iter = t)
  out
}

