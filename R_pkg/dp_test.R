# Testing for NaNs
library(glmgen)
set.seed(0)
n = 10000
y = exp(runif(n,0,185))

b = trendfilter(y,k=0,lambda=1e13)$beta
sum(b)

## # This should work: no NaNs (interruptions)
## for (r in 1:1000) {
##  y = exp(runif(n,0,185))
##  b = trendfilter(y,k=0,lambda=1e13)$beta
##  sum(b)
## } 

###
set.seed(100)
n = 10000
y = runif(n,0,1) # max(y)/min(y) ~500 for this draw

b = trendfilter(y/max(y),k=0,lambda=1e-35/max(y))$beta;
sum(b)

###
library(genlasso)
n = 10000
y = exp(runif(n,0,1))

# This should also "work" with the below specification,
# in the sense that no NaNs occur. But the errors to genlasso
# will just appear kind of large, since y is itself large
# y = exp(runif(n,0,185))

out = genlasso::fusedlasso1d(y)
b1 = out$beta
b2 = matrix(0,n,ncol(b1))
for (j in 1:ncol(b1)) {
  b2[,j] = glmgen::trendfilter(y,k=0,lambda=out$lambda[j])$beta
}
max(abs(b1-b2))

c1 = c2 = numeric(ncol(b1))
for (j in 1:ncol(b1)) {
  c1[j] = 0.5*sum((y-b1[,j])^2)+out$lambda[j]*sum(abs(diff(b1[,j])))
  c2[j] = 0.5*sum((y-b2[,j])^2)+out$lambda[j]*sum(abs(diff(b2[,j])))
}
sum(c1<c2)
min(c1-c2)
min(c2-c1)
