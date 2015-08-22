library(genlasso)
library(glmgen)
library(testthat)

y <- rnorm(100)

# test multiplication by "D"
for (k in 1:10)
  expect_equal(tfMultiply(y, k=k), as.numeric(getDtf(length(y),ord=k-1)%*%y))
