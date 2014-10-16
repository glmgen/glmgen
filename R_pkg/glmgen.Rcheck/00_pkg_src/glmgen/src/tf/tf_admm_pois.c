#include "tf.h"
#include "math.h"

static double b(double x)
{
  return exp(x);
}

static double b1(double x)
{
  return exp(x);
}

static double b2(double x)
{
  return exp(x);
}

void tf_admm_pois (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol)
{
  tf_admm_glm(y, x, n, k, max_iter, lam, beta, alpha, u, obj, iter,
     rho, obj_tol,
     &b, &b1, &b2);
}
