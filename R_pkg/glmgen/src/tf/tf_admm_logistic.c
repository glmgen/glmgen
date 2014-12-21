#include "tf.h"
#include "tf_glm_loss.h"

void tf_admm_logistic (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, int max_iter_admm, cs * DktDk)
{

  tf_admm_glm(y, x, w, n, k, max_iter, lam, beta, alpha, u, obj, iter,
     rho, obj_tol, max_iter_admm, DktDk, &logi_b, &logi_b1, &logi_b2);
}
