library(genlasso)
library(glmgen)
library(testthat)

set.seed(0)
EPS = 1e-4

n = 100
x = sort(runif(n, min=-2*pi, max=2*pi))
y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)

# Fused lasso w/o location values
outGenlasso = genlasso::trendfilter(y, ord=0)
outGlmgen = glmgen::trendfilter(y, k=0, lambda=outGenlasso$lambda)
expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)

# Fused lasso with location values
outGenlasso = genlasso::trendfilter(y, x, ord=0)
outGlmgen = glmgen::trendfilter(x, y, k=0, lambda=outGenlasso$lambda)
expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)

# Fused lasso using glmgen postions
outGenlasso = genlasso::trendfilter(y, x, ord=0)
outGlmgen = glmgen::trendfilter(x, y, k=0)
predGenlasso = coef(outGenlasso, lambda=outGlmgen$lambda)$beta
predGlmgen = predict(outGlmgen)
expect_true(abs(max(predGenlasso - predGlmgen)) < EPS)

p = trendfilter.control.list()
p$max_iter = 2000
p$obj_tol = 1e-10

# Higher order trendfiltering w/o location values
for (k in 1:2) {
  outGenlasso = genlasso::trendfilter(y, ord=k)
  outGlmgen = glmgen::trendfilter(y, k=k, lambda=outGenlasso$lambda, control=p)
  expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)
  print(abs(max(outGenlasso$beta - outGlmgen$beta)))
}

p$max_iter = 4000

# Higher order trendfiltering with location values
for (k in 1:2) {
  outGenlasso = genlasso::trendfilter(y, x, ord=k)
  outGlmgen = glmgen::trendfilter(x, y, k=k, lambda=outGenlasso$lambda, 
    control=p)
  expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < 1e-2)
  print(abs(max(outGenlasso$beta - outGlmgen$beta)))
}

# Polynomial signal test
lambda <- 1e4
n = 4
ks = c(0, 1, 2)
for (k in ks) {
  X = seq.int(0, n)
  y = X^k
  mod <- trendfilter(x = X, y = y, k = k, lambda = lambda)
  out <- predict(mod, lambda = lambda, type = "response")
  expect_true(abs(max(out - y)) < EPS)
  print(abs(max(out - y)))
}


