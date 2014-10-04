#include <stdio.h>

#include "cs.h"
#include "utils.h"
#include "tf.h"
#include "int_codes.h"

int main()
{

  int i;

  double * y;
  double * x;
  double * w;
  int n;
  int k;
  int family;
  int max_iter;
  int lam_flag;
  int obj_flag;
  double * lambda;
  int nlambda;
  double lambda_min_ratio;
  double * beta;
  double * obj;
  int * iter;
  double rho;
  double obj_tol;

  /* Input parameters; will usually be given by parent function to tf_admm */
  n = 8;
  k = 3;
  family = FAMILY_GAUSSIAN;
  max_iter = 25;
  lam_flag = 0;
  obj_flag = 1;
  nlambda = 3;
  lambda_min_ratio = 1e-4;
  rho = 1;
  obj_tol = 1e-12;

  y = (double *) malloc(n * sizeof(double));
  x = (double *) malloc(n * sizeof(double));
  w = (double *) malloc(n * sizeof(double));

  lambda = (double *) malloc(nlambda * sizeof(double));
  beta = (double *) malloc(n * nlambda * sizeof(double));
  obj = (double *) malloc(max_iter * nlambda * sizeof(double));
  iter = (int *) malloc(max_iter * nlambda * sizeof(int));

  for (i = 0; i < n; i++) y[i] = i + (i > 2)*3;
  for (i = 0; i < n; i++) x[i] = i;
  for (i = 0; i < n; i++) w[i] = 1;

  /* Call the tf_admm function */
  tf_admm(y, x, w, n, k, family, max_iter, lam_flag, obj_flag,
          lambda, nlambda, lambda_min_ratio, beta, obj, iter, rho, obj_tol);

  printf("\n---------- lambda -------------------------------\n");
  for (i = 0; i < nlambda; i++) printf("%f\n", lambda[i]);

  printf("\n---------- beta_1 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i]);

  printf("\n---------- beta_2 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i + n]);

  printf("\n---------- beta_3 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i + n*2]);


  /* Free the allocated arrays */
  free(y);
  free(x);
  free(w);
  free(lambda);
  free(beta);
  free(obj);
  free(iter);
}
