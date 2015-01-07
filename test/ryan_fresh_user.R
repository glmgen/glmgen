library(glmgen)

n = 1000
set.seed(0)
x = runif(n, min=-2*pi, max=2*pi)
y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
out = trendfilter(y, x, k=2)

# * why does is take so many iterations to converge at lambda_max??

matplot(1:nrow(out@obj),out@obj,type="l",lty=1)

xx = seq(-2*pi,2*pi,length=1000)
lambda = out@lambda[25]
f = predict(out,lambda=lambda)
ff = predict(out,x.new=xx,lambda=lambda)

plot(x,y)
points(out@x,f,pch=19,col=4)
lines(xx,ff,lwd=2,col=2)

# * how to flag someone that "thinning" has been turned on? when i ran it for 100 points
# and k=3, i went to predict with no x argument and got back predictions at 93 points
# ---> I added a warning flag, but then commented it out.  The user is just going to
# have to be aware and use out@x and out@y instead of x and y

# * why does thinning not return evenly spaced points?
# n = 100
# set.seed(0)
# x = runif(n, min=-2*pi, max=2*pi)
# y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
# out = trendfilter(y, x, k=3)
# diff(out@x)
# ---> I think it's because some bins are not occupied at all!

# * Right now thinning exits if m==1 in Veeru's code. I.e., if thinning were to reduce
# everything to one point, we quit. Shouldn't we verbalize this failure?

#########
# Mixture of Gaussians for x's!!

n = 1000
set.seed(0)
x = c(rnorm(n/2,mean=-pi,sd=1),rnorm(n/2,mean=pi,sd=1))
y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
out = trendfilter(y, x, k=2)

xx = seq(-2*pi,2*pi,length=1000)
lambda = out@lambda[30]
f = predict(out,lambda=lambda)
ff = predict(out,x.new=xx,lambda=lambda)

plot(x,y)
points(out@x,f,pch=19,col=4)
lines(xx,ff,lwd=2,col=2)

out2 = smooth.spline(x,y)
lines(xx,predict(out2,xx)$y,lwd=2,col=3)

#########
# Try the different link functions 

n = 1000
set.seed(0)
x = c(rnorm(n/2,mean=-pi,sd=1),rnorm(n/2,mean=pi,sd=1))
f0 = 1.5*sin(x) + sin(2*x) 
p = 1/(1+exp(-f0))
y = rbinom(n,1,p)
out = trendfilter(y, x, k=2, family="logistic")

xx = seq(-2*pi,2*pi,length=1000)
lambda = out@lambda[30]
f = predict(out,lambda=lambda)
ff = predict(out,x.new=xx,lambda=lambda)

plot(x,y)
o = order(x)
lines(x[o],p[o],col=3)
points(out@x,1/(1+exp(-f)),pch=19,col=4)
lines(xx,1/(1+exp(-ff)),lwd=2,col=2)

###

n = 1000
set.seed(0)
x = c(rnorm(n/2,mean=-pi,sd=1),rnorm(n/2,mean=pi,sd=1))
f0 = 1.5*sin(x) + sin(2*x) 
mu = exp(f0)
y = rpois(n,mu)
out = trendfilter(y, x, k=2, family="poisson")

xx = seq(-2*pi,2*pi,length=1000)
lambda = out@lambda[30]
f = predict(out,lambda=lambda)
ff = predict(out,x.new=xx,lambda=lambda)

plot(x,y)
o = order(x)
lines(x[o],mu[o],col=3)
points(out@x,exp(f),pch=19,col=4)
lines(xx,exp(ff),lwd=2,col=2)

