#!/usr/bin/env Rscript
library(genlasso)
n = 6
k = 0

x = 0:(n-1)
y = x + (x > 2)*3

D = getDtfSparse(n,k+1)
Dt = t(D)
DtD = t(D) %*% D
DDt = D %*% t(D)
Dy = D %*% y
b = solve(DDt, Dy)
kern = 0.5 * DDt + Diagonal(n-k-2)

cat("\n---------- D_(k+1) -----------------------------\n");
D
cat("\n---------- D_(k+1)^T ---------------------------\n");
Dt
cat("\n---------- D_(k+1)^T * D_(k+1) -----------------\n");
DtD
cat("\n---------- D_(k+1) * D_(k+1)^t -----------------\n");
DDt
cat("\n---------- y (for comparison) ------------------\n");
matrix(y,ncol=1)
cat("\n---------- D * y -------------------------------\n");
matrix(Dy,ncol=1)
cat("\n---------- b solution from -> DDt b = Dy -------\n");
matrix(b,ncol=1)
cat("\n---------- 0.5 * D_(k+1) * D_(k+1)^t + I_(k) ---\n");
kern
