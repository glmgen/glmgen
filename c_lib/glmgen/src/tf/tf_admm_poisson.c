#include "tf.h"
#include "tf_glm_loss.h"

void tf_admm_poisson (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter, 
       double rho, double obj_tol, cs * DktDk)
{
  tf_admm_glm(y, x, w, n, k, max_iter, lam, beta, alpha, u, obj, iter,
       rho, obj_tol, DktDk, &pois_b, &pois_b1, &pois_b2);
}
