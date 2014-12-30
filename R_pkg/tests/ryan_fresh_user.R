library(glmgen)

n = 1000
set.seed(0)
x = runif(n, min=-2*pi, max=2*pi)
y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
out = trendfilter(y, x, k=2)

# changed trendfiltering to trend filtering everywhere
# put in Veeru's full name everywhere, put my name last everywhere

# * should make the default k=2 in examples
# * why was the example trend filtering call OK with the different sized y and x??
# * what's up with out@p, why not p=1?
# * why have an objective flag in trendfilter?
# * why does is take so many iterations to converge at lambda_max??
# * how to know how many iters the alg took for each lambda?

matplot(1:nrow(out@obj),out@obj,type="l",lty=1)

xx = seq(-2*pi,2*pi,length=1000)
lambda = out@lambda[25]
f = predict(out,lambda=lambda)
ff = predict(out,x.new=xx,lambda=lambda)

plot(x,y)
#points(out@x,f,pch=19)
lines(xx,ff,lwd=2,col=2)

# * see also on trendfilter documentation is bogus
# * i think we shouldn't name it predict --- how to make see the help file?
# * should remove the error when predict is outside observed range

# * also, how to flag someone that "thinning" has been turned on? when i ran it for 100 points
# and k=3, i went to predict with no x argument and got back predictions at 93 points
# ---> I added a warning flag

# * why does thinning not return evenly spaced points?
# n = 100
# set.seed(0)
# x = runif(n, min=-2*pi, max=2*pi)
# y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
# out = trendfilter(y, x, k=3)
# diff(out@x)
# ---> I think it's because some bins are not occupied at all

# * I changed the entry point to thinning as follows:
# if (thinning && cond > control$x_cond)

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
#points(out@x,f,pch=19)
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
out = trendfilter(y, x, k=2, family="logistic", verb=T)

# i'd like to separate out max_iter for ADMM in the Gaussian case, and
# max_outer_iter for proximal Newton in the logistic & Poisson cases.
# This is because max_iter=200 is appropriate for the Gaussian case, but
# we're only going to want max_outer_iter=50 for the logistic & Poisson cases
# (and max_inner_iter=200). This is because the outer iterations are on a very
# different scale

xx = seq(-2*pi,2*pi,length=1000)
lambda = out@lambda[30]
f = predict(out,lambda=lambda)
ff = predict(out,x.new=xx,lambda=lambda)

plot(x,y)
o = order(x)
lines(x[o],p[o],col=3)
points(out@x,1/(1+exp(-f)),pch=19)
lines(xx,1/(1+exp(-ff)),lwd=2,col=2)

###
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
points(out@x,exp(f),pch=19)
lines(xx,exp(ff),lwd=2,col=2)

# Where is this message from? 
# "NaN for lambda[0], n=852, k=2"
# Found it; it means failure to converge at the first lambda value, we need to do
# better here.  Can we also return the status vector
