#include <math.h>

#include "tf.h"
#include "utils.h"

void tf_admm (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj, int * numiter,
              double rho, double eabs, double erel)
{
  beta[0] = 1;
  beta[1] = rho;
}

void tf_prime_dual (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj, int * numiter)
{
  beta[0] = 2;
}
