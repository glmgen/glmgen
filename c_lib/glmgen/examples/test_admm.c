#include "cs.h"
#include "utils.h"
#include "tf.h"
#include "test_utils.h"
#include "int_codes.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void test_admm_gauss();

int main()
{
  test_admm_gauss();
  return 0;
}
void test_admm_gauss(int mode)
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
  double * pred;
  double * obj;
  int * iter;  
  double rho;
  double obj_tol;
  double predict_zero_tol;

  /* Input parameters; will usually be given by parent function to tf_admm */
  n = 8;
  k = 3;
  family = FAMILY_GAUSSIAN;
  max_iter = 100;
  lam_flag = 1;
  obj_flag = 1;
  nlambda = 1;
  lambda_min_ratio = 1e-4;
  rho = 1;
  obj_tol = 1e-12;  

  y = (double *) malloc(n * sizeof(double));
  x = (double *) malloc(n * sizeof(double));
  w = (double *) malloc(n * sizeof(double));

  lambda = (double *) malloc(nlambda * sizeof(double));
  beta = (double *) malloc(n * nlambda * sizeof(double));
  pred = (double *) malloc(n * sizeof(double));
  obj = (double *) malloc(max_iter * nlambda * sizeof(double));
  iter = (int *) malloc( nlambda * sizeof(int));

  srand(5490);
  srand(time(NULL));

  x[0] = 0;
  for (i = 1; i < n; i++) x[i] = x[i-1] + ((rand() % 100)+1)/100.;
  for (i = 0; i < n; i++) y[i] = x[i] * x[i] + 0.5 * ((rand() % 100)+1)/100.;
  for (i = 0; i < n; i++) w[i] = 0;

  lambda[0] = 1;
  /* Call the tf_admm function */
  tf_admm(y, x, w, n, k, family, max_iter, lam_flag, obj_flag,
          lambda, nlambda, lambda_min_ratio, beta, obj, iter, rho, obj_tol);

  printf("\n--------------- (x,y) ---------------------------\n");
  for (i = 0; i < n; i++) printf("%.3f\t%1.f\n", x[i], y[i]);

  printf("\n---------- lambda -------------------------------\n");
  for (i = 0; i < nlambda; i++) printf("%f\n", lambda[i]);

  printf("\n---------- beta_1 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i]);

  /*
  printf("\n---------- beta_2 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i + n]);

  printf("\n---------- beta_3 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i + n*2]);
  */
  predict_zero_tol = 1e-12;
  
  double err;
  for(i = 0; i < nlambda; i++)
  {
    int offset = i*n;
    tf_predict_gauss(beta + offset, x, n, k, x, n, pred, predict_zero_tol);
    
    err = max_diff(beta + offset, pred, n);    
    printf("Prediction difference at input points=%E\n", err);
    if(!(err < 1e-12 ))
      printf("Prediction failed at input points (n=%d,k=%d,lam=%.4f,err=%E)\n",n,k,lambda[i], err);
  }
  /* Logistic loss */
  
  
  family = FAMILY_LOGISTIC;
  double bernouli_p, uniform;
  for (i = 0; i < n; i++)
  {
    bernouli_p = 1/ ( 1 + exp(-(x[i]-n/4.)*(x[i]-n/4.)) );
    uniform = (rand() % 101)/100.;    
    y[i] = uniform <= bernouli_p ? 1 : 0;
  }

  printf("\n--------------- (x,y) ---------------------------\n");

  for (i = 0; i < n; i++) printf("%.3f\t%1.f\n", x[i], y[i]);
  
  lambda[0] = 1; 
  tf_admm(y, x, w, n, k, family, max_iter, lam_flag, obj_flag,
          lambda, nlambda, lambda_min_ratio, beta, obj, iter, rho, obj_tol);
  
  printf("\n---------- lambda -------------------------------\n");
  for (i = 0; i < nlambda; i++) printf("%f\n", lambda[i]);
          
  printf("\n---------- beta_1 (logistic) -----------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i]);
  
  /* Free the allocated arrays */
  free(y);
  free(x);
  free(w);
  free(lambda);
  free(beta);
  free(obj);
  free(iter);
  free(pred);
  
}
