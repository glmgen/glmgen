#include "tf.h"
#include "math.h"

/* b(x) = log( 1 + exp(x) ).
 * Avoid computing exp(x) when x >> 0*/
static double b(double x)
{
  return x <= 0 ? log( 1 + exp(x) ) : x + log(1+exp(-x));
}

/* b1(x) = b'(x), the first derivative */
static double b1(double x)
{
  return x > 0 ? 1 / (1 + exp(-x) ) : exp(x)/ (1 + exp(x) );
  /* return 1. / (1 + exp(-x) ); */
}

/* b2(x) = b''(x), the second derivative */
static double b2(double x)
{
  x = -fabs(x);
  return exp(x-2*log(1+exp(x)));
  /* return exp(x) / ((1 + exp(x))*(1 + exp(x))); */
}

void tf_admm_logistic (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol)
{

  tf_admm_glm(y, x, n, k, max_iter, lam, beta, alpha, u, obj, iter,
     rho, obj_tol,
     &b, &b1, &b2);

}
