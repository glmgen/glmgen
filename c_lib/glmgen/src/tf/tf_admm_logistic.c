#include "tf.h"
#include "math.h"

static double b1(double x)
{
  return 1. / (1 + exp(-x) );
}

static double b2(double x)
{
  return exp(x) / ((1 + exp(x))*(1 + exp(x)));
}

void tf_admm_logistic (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR)
{
  
  tf_admm_glm(y, x, n, k, max_iter, lam, beta, alpha, u, obj,
    rho, obj_tol, sparseQR,
    &b1, &b2);
    
}
