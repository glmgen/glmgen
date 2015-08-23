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

# Higher order trendfiltering w/o location values
for (k in 1:2) {
  outGenlasso = genlasso::trendfilter(y, ord=k)
  outGlmgen = glmgen::trendfilter(y, k=k, lambda=outGenlasso$lambda)
  # expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)
  print(abs(max(outGenlasso$beta - outGlmgen$beta)))
}

# Higher order trendfiltering with location values
for (k in 1:2) {
  outGenlasso = genlasso::trendfilter(y, x, ord=k)
  outGlmgen = glmgen::trendfilter(x, y, k=k, lambda=outGenlasso$lambda)
  # expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)
  print(abs(max(outGenlasso$beta - outGlmgen$beta)))
}


