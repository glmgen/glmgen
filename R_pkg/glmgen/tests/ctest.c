#include "cs.h"
#include "utils.h"
#include "tf.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define PI 3.14159265

void test_admm_gauss(int n, int k);

double max(double a, double b)
{
  return a >= b ? a : b; 
}

double max_diff(double* a, double* b, int n)
{
  if(n<=0)
    return 0.0;
  double d = 0.0;
  while(n--)
  {
    d = max(d, fabs(a[n] - b[n]));
  }
  return d;
}

int main()
{
  int n=5;
  int k=0;
  test_admm_gauss(n, k);
  /*for(k = 1; k < 2; k++)
    {
    test_admm_gauss(n, k);
    }*/
  return 0;
}


void test_admm_gauss(int n, int k)
{
  int i;

  double * y;
  double * x;
  double * w;
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
  int * df;
  int * iter;
  int * status;
  double rho;
  double obj_tol;
  double predict_zero_tol;
  int verb;
  int tridiag;

  verb = 1;
  family = FAMILY_GAUSSIAN;
  tridiag = 0;
  max_iter = 100; // admm iters
  lam_flag = 1;
  obj_flag = 1;
  nlambda = 1;
  lambda_min_ratio = 1e-4;
  rho = 1;
  obj_tol = 1e-8;
  double alpha_ls = 0.5, gamma_ls = 0.9;
  int max_iter_ls = 10, max_iter_newton = 1;
  int verbose = 1;

  y = (double *) malloc(n * sizeof(double));
  x = (double *) malloc(n * sizeof(double));
  w = (double *) malloc(n * sizeof(double));

  lambda = (double *) malloc(nlambda * sizeof(double));	
  df = (int *) malloc(nlambda * sizeof(int));
  beta = (double *) malloc(n * nlambda * sizeof(double));
  pred = (double *) malloc(n * sizeof(double));   
  iter = (int *) malloc( nlambda * sizeof(int));
  status = (int *) malloc( nlambda * sizeof(int));
  int max_iter_outer = (family == FAMILY_GAUSSIAN) ? max_iter : max_iter_newton;
  obj = (double *) malloc(max_iter_outer * nlambda * sizeof(double));

  lambda[0] = 1e4;
  for (i = 0; i < n; i++) x[i] = i;
  for (i = 0; i < n; i++) y[i] = pow(x[i], (double) k);
  for (i = 0; i < n; i++) w[i] = 1.0;

  double* beta0 = NULL;
  /* Call the tf_admm function */
  tf_admm(x,y,w,n,k,family,
      max_iter,beta0,lam_flag,lambda,
      nlambda,lambda_min_ratio, tridiag, df,
      beta,	obj, iter, status,
      rho, obj_tol, obj_tol, alpha_ls, gamma_ls,
      max_iter_ls, max_iter_newton, verbose);

  /*  printf("\n--------------- (x,y) ---------------------------\n");*/
  /*  for (i = 0; i < n; i++) printf("%.2f\t%.2f\n", x[i], y[i]);*/

  /*  printf("\n---------- lambda -------------------------------\n");*/
  /*    for (i = 0; i < nlambda; i++) printf("%f\n", lambda[i]);*/

  /*  printf("\n---------- beta_1 -------------------------------\n");*/
  /*    for (i = 0; i < n; i++) printf("%f\n", beta[i]);*/


  /* Prediction */
  
  predict_zero_tol = 1e-12;

  double err;
  for(i = 0; i < nlambda; i++)
  {
  int offset = i*n;
  tf_predict_gauss(beta + offset, x, n, k, x, n, pred, predict_zero_tol);

  if(verb)
  {
    printf("predicted values --------\n");
    print_array(pred, n);
    printf("expected values --------\n");
    print_array(beta + offset, n);
  }
  err = max_diff(beta + offset, pred, n);
  printf("Prediction difference at input points=%E\n", err);
  if(!(err < 1e-8 ))
  printf("Prediction failed at input points (n=%d,k=%d,lam=%.4f,err=%E)\n",n,k,lambda[i], err);
}
     

  /* Logistic loss */
  /*

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
     lambda, nlambda, lambda_min_ratio, beta, obj, iter, status, rho, obj_tol);

     printf("\n---------- lambda -------------------------------\n");
     for (i = 0; i < nlambda; i++) printf("%f\n", lambda[i]);

     printf("\n---------- beta_1 (logistic) -----------------------\n");
     for (i = 0; i < n; i++) printf("%f\n", beta[i]);
     */
  /* Free the allocated arrays */
  free(y);
  free(x);
  free(w);
  free(lambda);
  free(beta);
  free(obj);
  free(df);
  free(iter);
  free(status);
  free(pred);

}



// cd .. && sudo R CMD INSTALL glmgen && cd - && \
// cp src/glmgen.so src/libglmgen.so && \
// gcc -g tests/ctest.c -I inst/include/ -Lsrc/ -lm -lglmgen && gdb a.out
