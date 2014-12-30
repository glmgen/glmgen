
library(R.matlab)

cvx_data = readMat("../../glmcvx/crit_logi_cvx.mat");
crit_cvx = cvx_data$criterions

load("cvx_cmp.RData"); # from prox_newton

pdf(file="cvx_cmp.pdf", height=4, width=4*length(kvals))
for(n in nvals) {
  par(mfrow=c(1,length(kvals)))

  for(k in kvals) {
    cat("k=",k,", n=", n, "\n")
    a = crit_a[!is.na(match( crit_a[,1],n) & match( crit_a[,2], k)),]
    c = crit_cvx[!is.na(match( crit_cvx[,1],n) & match( crit_cvx[,2], k)),]

    plot(1:nrow(a), a[,5], type="l", xlab="lam num", ylab="optimal value")
    lines(1:nrow(a), c[,5], col=2)
    legend("topright",legend=c("ADMM","CVX"),
           col=c("black","red"),pch=21)

    title(paste("n=", n, ", k=", k))

  }
}
dev.off()
