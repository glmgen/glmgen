library(glmgen)

n = 10000
lam = 1e13
y = exp(runif(n,0,185))
sum(is.nan(trendfilter(y,k=0,lambda=lam)$beta))
sum(is.nan(trendfilter(y/max(y),k=0,lambda=1e13)$beta))
sum(is.nan(trendfilter(y/max(y),k=0,lambda=1)$beta))
