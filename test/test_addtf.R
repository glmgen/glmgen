
library(glmgen)
source("funs.R")
source("addtf.R")

# Setup the roblem
set.seed(2)
n = 1000;
d = 2;

# Spaces out x to avoid thinning
space = function(x) {
  out = sort(x, index.return=T)
  x = out$x; 
  ord = order( out$ix )

  n = length(x)
  dx = diff(x)
  delta = 1/n
  for(i in 1:(n-1)) {
    x[i+1] = x[i] + dx[i] +  delta
  }

  x[ord]
}

x = matrix( runif(n*d,0,1), ncol=d)
for(j in 1:d) {
  x[,j] = space( x[,j] )
}

f = x
f[,1] = 5*sin( x[,1] * pi/2 )
f[,2] = exp( (x[,2] ) )
y = rowSums(f) + 10 + rnorm(n,sd=0.1)
lambda = 1
k = 1
ks = rep(k,d)

# Backfitting settings
max_iter = 10

# Backfitting
out = addtf(y=y, x=x, ks=ks, lambda=lambda, verbose=T)

# Plot
beta = out$beta
obj = out$obj
pdf("~/Dropbox/addmod/plots/sin_exp_lam_1.pdf", height=4, width=4)

par(mfrow=c(1,1))
for(j in 1:d) {
  ord = order(x[,j])
  plot(x[ord,j],f[ord,j],type="l",col="blue",
       xlab=paste("x[",j,"]"),
       ylab="f/fit",ylim=range( c(f[,j],beta[,j]), na.rm=TRUE))
  lines(x[ord,j],beta[ord,j],col="red")
  legend("topright", legend=c("f", "fit"), col = c("blue", "red"), pch="-")

  title(paste("n=", n, ", k=", k, ", j=", j, sprintf(",lam=%g", lambda)))
}

plot(1:max_iter, obj, xlab="#iter", ylab="objective", ylim = range(obj, na.rm=TRUE))

dev.off()

