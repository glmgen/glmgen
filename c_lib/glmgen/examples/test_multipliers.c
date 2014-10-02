#include "cs.h"
#include "utils.h"
#include "tf.h"
#include <string.h>
#include <tgmath.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>

void test_mult();
void test_pred();

int main() 
{
  
  test_mult();
  test_pred();
  
  return 0;
}

int test_d(int n, int k, double *x, double *a, int verbose);
double max_diff(double *x, double *y, int n);
void test_mult()
{
  double *x, *a;

  int n, verb, rep, reps, num_failed;
  n = 8;  
  x = (double*)malloc(n*sizeof(double));
  a = (double*)malloc(n*sizeof(double));

  verb = 0;
  reps = 10;
  num_failed = 0;  
  
  srand(5489);

  for(rep=0; rep<reps; ++rep) {
    int i=0;
    x[0] = 0;
    for(i=0; i<n;++i) {
      a[i] = i*i*i + (rand() % 10) * i;
      if(i)
        x[i] = x[i-1] + (rand() % 100) / 100.;
    }
    num_failed += test_d(n, 1, x, a, verb);
    num_failed += test_d(n, 2, x, a, verb);
    num_failed += test_d(n, 3, x, a, verb);
    num_failed += test_d(n, 7, x, a, verb);
  }

  printf("\nNumber of failed tests: %d/%d\n", num_failed, reps * 4 );

  free(x);
  free(a);
}

void test_pred()
{
  double *x, *x0, *beta;

  int n, n0, verb, rep, reps, num_failed;

  n = 8;
  n0 = 5;

  if(n0 >= n)
    n0 = n-1;
  
  x = (double*)malloc(n*sizeof(double));
  x0 = (double*)malloc(n0*sizeof(double));
  beta = (double*)malloc(n*sizeof(double));
  
  verb = 1;
  reps = 10;
  num_failed = 0;
  
  srand(5489);
  
  for(rep=0; rep<reps; ++rep) {
    int i=0;
    x[0] = 0;
    for(i=1; i<n;++i) 
    {
      x[i] = x[i-1] + (rand() % 100) / 100.;
    }
    for(i=0; i<n0; i++)
    {
      x0[i] = (x[i] + x[i+1]) / 2;
    }
    
    for(i=0; i<n; i++)
    {
      beta[i] = i*i;
    }    
    
    num_failed += test_predict(x,n,2,x0,n0,beta,verb);
  }

  printf("\nNumber of failed tests: %d/%d\n", num_failed, reps);

  free(x);
  free(x0);
  free(beta);
}
void print_array(double *x, size_t n);


int test_d(int n, int k, double *x, double *a, int verb)
{
  double *b, *b1;
  cs *D, *Dt;

  b = (double*)calloc(n, sizeof(double));
  b1 = (double*)calloc(n-k,sizeof(double));

  printf("\nTesting n=%d, k=%d--------------\n\n", n, k);
  if(verb) 
  {
    printf("x-------------------------------\n");
    print_array(x,n);

    printf("a-------------------------------\n");
    print_array(a,n);
  }
  
  tf_dx(x,n,k,a,b);
  if(verb)
  {
    printf("D^k*a (in place)-----------------\n");
    print_array(b,n-k);
  }
  D = tf_calc_dk(n, k-1, x); /* D^k, not D^(k-1) */  
  if(verb)
  {
    printf("D^%d-------------------------\n", k);
    cs_print(D, 0); 
  }

  cs_gaxpy(D,a,b1);
  if(verb)
  {
    printf("D^k*a (multiplication)-----------\n");
    print_array(b1,n-k);
  }

  double err = max_diff(b, b1, n-k);
  printf("Dx err=%f\n", err);
  if(err > 1e-10 * pow(10.0,k/2)) 
  {
    printf("*************D Test failed(n=%d,k=%d,err=%E)**********\n",n,k,err);
    return 1;
  }

  b = (double*)realloc(b, n*sizeof(double));
  tf_dtx(x,n,k,a,b);
  if(verb)
  {
    printf("(D^k)'*a (in place)-----------------\n");
    print_array(b,n);
  }

  free(b1);
  b1 = (double*)calloc(n,sizeof(double));
  Dt = cs_transpose(D,1);
  cs_gaxpy(Dt,a,b1);
  if(verb)
  {
    printf("(D^k)'*a (multiplication)-----------\n");
    print_array(b1,n);
  }

  err = max_diff(b, b1, n);
  printf("Dtx err=%f\n", err);
  if(err > 1e-10* pow(10.0,k/2)) 
  {
    printf("*************Dt Test failed(n=%d,k=%d,err=%E)**********\n",n,k,err);
    return 1;
  }

  free(b);
  free(b1);
  cs_spfree(D);
  cs_spfree(Dt);
  
  return 0;
}
void predict_gauss_explicit(double * beta, double * x, int n, int k, 
	double * x0, int n0, double * pred,
	double zero_tol);
	
int test_predict(double *x, int n, int k, double *x0, int n0, double *beta, int verb)
{

  double *pred_exp = (double*)calloc(n0, sizeof(double));
  double *pred = (double*)calloc(n0, sizeof(double));
  
  double zero_tol = 1e-10;
  tf_predict_gauss(beta, x, n, k, x0, n0, pred, zero_tol);
  /* predict_gauss_explicit(beta, x, n, k, x0, n0, pred_exp, zero_tol); */
  
  double err = max_diff(pred, pred_exp, n0);
  
  printf("predict err=%f\n", err);
  if(err > 1e-12) 
  {
    printf("*************Predict Test failed(n=%d,k=%d,err=%E)**********\n",n,k,err);
    return 1;
  }
  
  free(pred_exp);
  free(pred);
  return 0;  
}
void print_array(double *x, size_t n)
{
  int i=0;

  for(i=0; i<n; i++)
  {
    printf("%f\n", x[i]);
  }
}

double norm(double * x, int n) {
  assert(x);
  double sum = 0.;
  int i=0;
  for(i=0; i<n; ++i)
    sum += (x[i] * x[i]);

  return sqrt(sum);
}
double residual(cs * A, double * x, double * b) {

  int m = A->m;
  double * neg_b = (double*)malloc(m*sizeof(double));
  int i=0;
  for(i=0; i < m; ++i)
    neg_b[i] = -b[i];
  cs_gaxpy(A,x,neg_b);

  double res = norm(neg_b,m);

  cs_free(neg_b);
  return res;
}

double max_diff(double *x, double *y, int n)
{
  int i=0;
  double e=0;
  for(i=0; i<n; i++)
    if( fabs(x[i] - y[i]) > e)
      e = fabs(x[i] - y[i]);

  return e;
}

static void poly_coefs(double *x, int n, int k, 
  double *beta, double *phi)
{
  int i;

  double *u = (double *)malloc((n)*sizeof(double));
  
  phi[0] = beta[0];
  
  for(i=1; i<k+1; i++)
  {
    tf_dx(x,n,i,beta,u);
    phi[i] = u[0] / ( (x[i]-x[0]) * glmgen_factorial(i-1));
  }
  
  free(u);
}

void predict_gauss_explicit(double * beta, double * x, int n, int k, 
	double * x0, int n0, double * pred,
	double zero_tol)
{
  int i=0, j=0;

  double *phi = (double *)malloc((k+1)*sizeof(double));
  poly_coefs(x,n,k,beta,phi);

  double *theta = (double *)malloc((n)*sizeof(double));
  tf_dx(x,n,k+1,beta,theta);

  for (i=0; i<n-k-1; i++) if (fabs(theta[i])<zero_tol) theta[i]=0;
  
  /* Compute the predictions at each new point x0 */
  double h;
  for (j=0; j<n0; j++) {
    pred[j] = 0;
   
    /* Loop over x points, polynomial basis */
    for (i=0; i<k+1; i++) {  
      h = 1;
      int l=0;
      for (l=0; l<i; l++) {
        h *= (x0[j]-x[l]);
      }
      pred[j] += phi[i]*h;
    }

    /* Loop over x points, falling fact basis */
    for (i=0; i<n-k-1; i++) {
      /* If the current x0 is too small, then break */
      if (x0[j]<x[i+k]) break;

      if (theta[i]!=0) {
        h = 1;
        int l=0;
        for (l=0; l<k; l++) {
          h *= (x0[j]-x[i+l+1]);
        }
        pred[j] += theta[i]*h;
      }
    }
  }

  free(phi);
  free(theta);    
}
