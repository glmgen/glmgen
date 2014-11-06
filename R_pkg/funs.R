library(Matrix)

Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

getDtfSparse = function(n,ord) {
  D = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
  D0 = D
  for (i in Seq(1,ord)) D = D0 %*% D
  return(D[Seq(1,n-ord-1),]) 
}

getDtf = function(n, ord) {
  return(as.matrix(getDtfSparse(n,ord)))
}

getT = function(n,k,x=1:n) {
  T = x
  if (k==0) T = T[-1]
  else if (k%%2==1) T = T[-c(1:((k+1)/2),(n-(k-1)/2):n)]  
  else if (k%%2==0) T = T[-c(1:((k+2)/2),(n-(k-2)/2):n)]
  return(T)
}

getG = function(n,k,x=1:n) {
  T = getT(n,k,x)
  G = matrix(0,n,n)
  G[,1] = rep(1,n)
  for (j in 2:n) {
    if (j<=k+1) {
      G[,j] = x^(j-1)
    }
    else {
      G[,j] = pmax(x-T[j-k-1],0)^k
    }
  }
  return(G)
}

getH = function(n,k,x=1:n) {
  if (k==0) return(lower.tri(matrix(0,n,n),diag=TRUE)+0)
  H = matrix(0,n,n)
  H[,1] = rep(1,n)
  for (j in 2:n) {
    if (j<=k+1) {
      H[,j] = apply(matrix(rep(x,each=j-1),ncol=n)-
         matrix(rep(x[1:(j-1)],n),ncol=n),2,prod)/factorial(j-1)
    }
    else {
      H[,j] = apply(pmax(matrix(rep(x,each=k),ncol=n)-
         matrix(rep(x[j-k-1+1:k],n),ncol=n),0),2,prod)/factorial(k)
    }
  }
  return(H)
}

getHslow = function(n,k) {
  D = matrix(0,n,n)
  D[1,1] = 1
  for (i in 2:(k+1)) D[i,] = getDtf(n,i-2)[1,]
  D[(k+2):n,] = getDtf(n,k)
  return(round(solve(D),3))
}

getHscaled = function(n,k,x=1:n) {
  H = getH(n,k,x)
  return(H*factorial(k))
}

compare = function(a,aa,lam) {
  f = predict(a,lambda=lam)$fit
  if ("trendfilter" %in% class(aa)) ff = coef(aa,lambda=lam)$beta
  else ff = predict(aa,lambda=lam)$fit
  n = length(a$y)
  plot(1:n,a$y)
  lines(1:n,f,col="blue")
  lines(1:n,ff,col="red")
  legend("bottomleft",col=c("blue","red"),
         legend=c("LARS","TF"),lty=1)
  return(max(abs(f-ff)))
}

kspline = function(k=0, n=100, seed=0, numknots=5, knots=NULL, weights=NULL) {
  if (is.null(knots)) knots=sample(1:n,numknots)
  else numknots = length(knots)

  if (is.null(weights))
    weights=sample(c(-1,1),numknots+k+1,replace=TRUE)*rnorm(numknots+k+1,2,1)
  
  beta = rep(0,n)
  for (i in 0:k) {
    beta = beta + weights[i+1]*(1:n)^k/n^k
  }
  for (j in 1:numknots) {
    x = 1:n-knots[j]
    x = x^k*(x>0)
    beta = beta + weights[j+k+1]*x/max(x)
  }

  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  y = beta + rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

pieconst = function(n=100, seed=0, numknots=5) {
  set.seed(seed)
 
  d = floor(n/numknots)
  beta = rep(sample(1:10,numknots),
    times=c(rep(d,numknots-1),n-(numknots-1)*d))
  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  
  y = beta + rnorm(n,sd=1)
  return(list(beta=beta,y=y))
}

pielin = function(n=100, seed=0) {
  set.seed(seed)
  n = 100
  knots = matrix(0,6,2)
  knots[,1] = c(1,20,45,60,85,100)
  knots[,2] = c(1,2,6,8,5,6)
  beta = rep(0,100)
  for (i in 1:(nrow(knots)-1)) {
    for (j in knots[i,1]:(knots[i+1,1])) {
      beta[j] = (knots[i+1,1]-j)/(knots[i+1,1]-knots[i,1])*knots[i,2] +
        (j-knots[i,1])/(knots[i+1,1]-knots[i,1])*knots[i+1,2]
    }
  }

  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  y = beta+rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

pielinx = function(x) {
  knots = matrix(0,6,2)
  knots[,1] = c(0,20,45,60,85,100)/100
  knots[,2] = c(1,2,6,8,5,6)
  ii = max(which(x>=knots[,1]))
  if (ii==nrow(knots)) {
    return(knots[nrow(knots),2])
  }
  else {
    a = (x-knots[ii,1])/(knots[ii+1,1]-knots[ii,1])
    return(knots[ii,2]*(1-a) + knots[ii+1,2]*a)
  }
}


piequad = function(n=100, seed=0) {
  set.seed(seed)
  knots = matrix(0,4,2)
  knots[,1] = c(1,33,60,100)
  knots[,2] = c(8,6,2,4)
  beta = rep(0,n)
  endval = 0
  for (i in 1:(nrow(knots)-1)) {
    sgn = (-1)^(i+1)
    mid = (knots[i,1]+knots[i+1,1])/2
    dif = knots[i+1,1]-knots[i,1]
    
    j = knots[i,1]
    intcp = endval - (sgn*n/5*(j-mid)^2/dif^2 + (knots[i+1,1]-j)/dif*knots[i,2] +
      (j-knots[i,1])/dif*knots[i+1,2])
    
    for (j in knots[i,1]:(knots[i+1,1])) {
      beta[j] = intcp + sgn*n/5*(j-mid)^2/dif^2 + (knots[i+1,1]-j)/dif*knots[i,2] +
        (j-knots[i,1])/dif*knots[i+1,2]
    }
    endval = beta[j]
  }
  
  beta = 10*(beta-min(beta))/(max(beta)-min(beta))
  y = beta + rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

piecub = function(n=100, seed=0) {
  set.seed(seed)
  n=100
  beta=rep(0,100)
  beta[1:40]=(1:40-20)^3
  beta[40:50]=-60*(40:50-50)^2 + 60*100+20^3
  beta[50:70]=-20*(50:70-50)^2 + 60*100+20^3
  beta[70:100]=-1/6*(70:100-110)^3 + -1/6*40^3 + 6000
  beta=-beta
  
  beta=10*(beta-min(beta))/(max(beta)-min(beta))
  y=beta+rnorm(n,0,1)

  return(list(beta=beta,y=y))
}

piecubx = function(x) {
  x = x*100
  if (x<=40) return((x-20)^3)
  if (x<=50) return(-60*(x-50)^2 + 60*100+20^3)
  if (x<=70) return(-30*(x-50)^2 + 60*100+20^3)
  if (x<=100) return(-1/6*(x-110)^3 + -1/6*40^3 + 6000)
}

smoothwiggly.fun = function(x) {
  f = function(a,b,c,x) return(a*(x-b)^2+c)
  fp = function(a,b,c,x) return(2*a*(x-b))

  a=-1; b=1/4; c=1;  
  if (x<=1/3) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=1/3;
  a=1; b=xx-fp(aa,bb,cc,xx)/(2*a); c=f(aa,bb,cc,xx)-a*(xx-b)^2;  
  if (x<=2/3) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=2/3;
  b=0.7; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.775) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.775;
  b=0.8; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.825) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.825;
  b=0.85; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.875) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.875;
  b=0.9; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.925) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.925;
  b=0.95; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.975) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.975;
  b=1; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  return(f(a,b,c,x))
}

smoothwiggly = function(n=100, x=1:n/n) {
  u = rep(0,n)
  for (i in 1:n) u[i]=smoothwiggly.fun(x[i])
  return(u)
}

maxlam2 = function(y,k) {
  n = length(y)
  D = getDtfSparse(n,k)
  x = qr(t(D))
  u = backsolveSparse(x,y)
  return(list(u=u,maxlam=max(abs(u))))
}

## maxlam3 = function(y,k) {
##   n = length(y)
##   D = getDtfSparse(n,k)
##   x = qr(crossprod(t(D)))
##   u = backsolveSparse(x,D%*%y)
##   return(list(u=u,maxlam=max(abs(u))))
## }

backsolveSparse = function(QR, b) {
  #R = qrR(QR)
  #x = solve(R, qr.qty(QR,b)[Seq(1,nrow(R))])
  #return(as.numeric(x))
  R = qr.R(QR)
  x = solve(R, qr.qty(QR,b)[Seq(1,nrow(R))])
  if (length(QR@q)==0) return(x)
  else return(x[Order(QR@q+1)])
}

Order = function(x) {
  n = length(x)
  o = numeric(n)
  o[x] = Seq(1,n)
  return(o)
}

crit = function(y,k,lambda,beta) {
  return(0.5*sum((y-beta)^2) + lambda*sum(abs(diff(beta,differences=k+1))))
}

bx = function(x,lam) {
  y = x
  y[x>lam] = lam
  y[x< -lam] = -lam
  return(y)
}

st = function(x,s,r) {
  z = rep(0,length(x))
  z[x>s] = x[x>s]-s
  z[x< -s] = x[x< -s]+s
  # leave the first r coordinates untouched
  z[1:r] = x[1:r]
  return(z)
}

## Objective
objective = function(x, y, k, lambda, beta) {
  Db = diff(beta);
  if( k >= 1 ) {
    for( i in 1:k) {
      Db = diff( i * (1/diff(x,i)) * Db)
    }
  }

  0.5 * sum((y-beta)*(y-beta)) + lambda * sum(abs(Db))
}

objectives = function(x, y, k, lam, beta ) {
  nlam = length(lam);
  obj = numeric(nlam);
  for(i in 1:nlam) {
    obj[i] = objective(x, y, k, lam[i], beta[,i]); 
  }
  obj
}

critD = function(y, D, lam, beta) {
  return (sum((y-beta)^2)/2 + lam*sum(abs(D %*% beta)))
}


