library("genlasso")
library("glmgen")
source("funs.R")

objective = function(x, y, k, lambda, beta) {
  Db = diff(beta);
  if( k > 1 ) {
    for( i in 1:k) {
      Db = diff( i * (1/diff(x,i)) * Db)
    }
  }

  0.5 * sum((y-beta)*(y-beta)) + lambda * sum(abs(Db))
}

kvals = c(1,2,3)
nvals = c(50,1e3,1e5)
xpows = c(1, 1/4)

# Test glmgen versus genlasso; uses varying x inputs
set.seed(42)
output01 = matrix(NA, ncol=6, nrow=0)
for (k in kvals) {
  for (n in nvals) {
    for (p in xpows) {
      x = (seq(0,1,length.out=n))^p
      y = sin(x * (2*pi)) + rnorm(n, sd=0.2)

      path = genlasso:::trendfilter(y, x, ord=k, maxsteps=20L)
      admm = glmgen::trendfilter(y, x, k=k, lambda = path$lambda,
                                          maxiter=200L, control=list(obj_tol=0))

      outtemp = cbind(k, n, p, 1:20L, rep(NA, 20L), rep(NA, 20L))
      colnames(outtemp)[4:6] = c("step", "path", "admm")
      for (i in 1:20L) {
        outtemp[i,5] = objective(x, y, k, path$lambda[i], path$beta[,i])
        outtemp[i,6] = objective(x, y, k, path$lambda[i], admm$beta[,i])
      }
      output01 = rbind(output01, outtemp)
    }
  }
}

# Test both admm algorithms; uses fixed input x
set.seed(42)
output02 = matrix(NA, ncol=6, nrow=0) # objective
output03 = matrix(NA, ncol=6, nrow=0) # timing
colnames(output02) = c("k", "n", "p", "step", "old", "new")
colnames(output03) = c("k", "n", "p", "step", "old", "new")
for (k in kvals) {
  for (n in nvals) {
    x = (seq(0,n-1,length.out=n))
    y = sin(x * (2*pi)) + rnorm(n, sd=0.2)

    t0 = system.time({
      admm_old = trendfilterSim::trendfilter(y=y, k=k, eabs=0, erel=0, maxiter=200L)
    })
    t1 = system.time({
      admm_new = glmgen::trendfilter(y, x, k=k, maxiter=200L, lambda=admm_old$lambda, control=list(obj_tol=0))
    })

    output03 = rbind(output03, c(k, n, 1, 1, t0[["elapsed"]], t1[["elapsed"]]))

    outtemp = cbind(k, n, 1, 1:20, rep(NA, 20), rep(NA, 20))
    colnames(outtemp)[4:6] = c("step", "old", "new")
    for (i in 1:20) {
      outtemp[i,5] = objective(x, y, k, admm_old$lambda[i], admm_old$beta[,i])
      outtemp[i,6] = objective(x, y, k, admm_old$lambda[i], admm_new$beta[,i])
    }
    output02 = rbind(output02, outtemp)
  }
}

output01[!is.na(match(output01[,4], c(1,10,20))),]
output02[!is.na(match(output02[,4], c(1,10,20))),]
output03


