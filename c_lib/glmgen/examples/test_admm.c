#include "cs.h"
#include "utils.h"
#include "tf.h"
#include "test_utils.h"
#include "int_codes.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define PI 3.14159265

void test_admm_gauss(int n, int k);

int main()
{
  int n=25;
  int k=3;
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
  double * alpha;
  double * z;
  double * obj;

  int * iter;
  int * status;
  double rho;
  int vary_rho;
  double obj_tol;
  double predict_zero_tol;
  int verb;

  verb = 0;
  /* Input parameters; will usually be given by parent function to tf_admm */
  family = FAMILY_GAUSSIAN;
  max_iter = 100;
  lam_flag = 0;
  obj_flag = 1;
  nlambda = 5;
  lambda_min_ratio = 1e-4;
  rho = 1;
  vary_rho = 1;
  obj_tol = 1e-12;

  y = (double *) malloc(n * sizeof(double));
  x = (double *) malloc(n * sizeof(double));
  w = (double *) malloc(n * sizeof(double));


  lambda = (double *) malloc(nlambda * sizeof(double));
  alpha = (double *) malloc(n * sizeof(double));
  z = (double *) malloc(n * sizeof(double));
  beta = (double *) malloc(n * nlambda * sizeof(double));
  pred = (double *) malloc(n * sizeof(double));
  obj = (double *) malloc(max_iter * nlambda * sizeof(double));
  iter = (int *) malloc( nlambda * sizeof(int));
  status = (int *) malloc( nlambda * sizeof(int));


  srand(5);
  /* srand(time(NULL)); */

  x[0] = 0;
  /* for (i = 1; i < n; i++) x[i] = x[i-1] + ((rand() % 100)+1)/100.; */
  for (i = 1; i < n; i++) x[i] = i;
  for (i = 0; i < n; i++) y[i] = sin(x[i] * 3.*PI/n);
 /* for (i = 0; i < n; i++) y[i] = x[i] * x[i] + 0.5 * ((rand() % 100)+1)/100.;*/
  for (i = 0; i < n; i++) w[i] = 1;
  
  for (i = 0; i < n; i++) z[i] = exp( rand() % 80 ) + exp( rand() % 45 );
/*  for (i = 0; i < n; i++) z[i] = sin(x[i] * 3.*PI/n);;*/
  for(i=0; i < n; i++) z[i] = z[i] * exp( -75 );
  tf_dp(n-k,z,1e-17, alpha);

  printf("alpha > 1e40 at\n");
  for(i = 0; i < n; i++) {
    if( fabs(alpha[i]) > 1e35 )
      printf("%d:%g  ", i, alpha[i]);
  }
  printf("\n");
  printf("z=%g, alpha=%g\n", l2norm(z,n-k), l2norm(alpha, n-k));

  for (i = 0; i < n; i++) printf("%g  ", z[i]);
  printf("\n");
  /* Call the tf_admm function */

  /* tf_admm(y, x, w, n, k, family, max_iter, lam_flag, obj_flag,
          lambda, nlambda, lambda_min_ratio, beta, obj, iter, status, rho, vary_rho, obj_tol);

  printf("\n--------------- (x,y) ---------------------------\n");
  for (i = 0; i < n; i++) printf("%.2f\t%.2f\n", x[i], y[i]);

  
  printf("\n---------- lambda -------------------------------\n");
  for (i = 0; i < nlambda; i++) printf("%f\n", lambda[i]);

  printf("\n---------- beta_1 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i]);

  printf("\n---------- beta_2 -------------------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", beta[i + n]);
*/

  /* Prediction */
  /*
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
    if(!(err < 1e-12 ))
      printf("Prediction failed at input points (n=%d,k=%d,lam=%.4f,err=%E)\n",n,k,lambda[i], err);
  }
  */

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
  free(iter);
  free(status);
  free(pred);

}
