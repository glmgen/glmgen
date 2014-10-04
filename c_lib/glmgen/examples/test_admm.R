library(genlasso)
n = 20
k = 2
nlambda = 3
lambda_min_ratio = 1e-4
rho = 1

# here, ord <-> k-1; hence my confusion in the C code:
D1 = getDtfSparse(n,ord=k)
D0 = getDtfSparse(n,ord=k-1)

x = 0:(n-1)
y = x + (x > 2)*3;

out = trendfilter(y = y, pos = x, ord = k)
max_lam = out$lambda[[1]]
min_lam = max_lam * lambda_min_ratio
lambda = exp((log(max_lam) * (nlambda - 0:(nlambda-1)) + log(min_lam) * 0:(nlambda-1)) / nlambda)

beta = predict(out,lambda=lambda)$fit

beta_max = beta[,1]
alpha_max = D0 %*% beta[,1]
u_max = solve( (D0 %*% t(D0)), D0 %*% (beta_max-y) / (rho * lambda[[1]]))

solve( (D0 %*% t(D0))) %*% D0 %*% (beta_max-y) / rho
# So, these should match from the c code:
{
  cat("---------- beta_max -------------------------------\n")
  cat(beta_max, sep="\n")
  cat("---------- alpha_max -------------------------------\n")
  cat(as.matrix(alpha_max), sep="\n")
  cat("---------- u_max -------------------------------\n")
  cat(as.matrix(u_max), sep="\n")

  cat("---------- lambda -------------------------------\n")
  cat(lambda, sep="\n")
  cat("---------- beta_1 -------------------------------\n")
  cat(beta[,1], sep="\n")
  cat("---------- beta_2 -------------------------------\n")
  cat(beta[,2], sep="\n")
  cat("---------- beta_3 -------------------------------\n")
  cat(beta[,3], sep="\n")
}