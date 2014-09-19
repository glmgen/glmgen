#include "cs.h"
#include "utils.h"
#include "tf.h"
#include <string.h>
#include <tgmath.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>

void test_d();
void print_array(double *x, size_t n);
double max_diff(double *x, double *y, int n);

int main() {

  double *x, *a;

  int n=8;
  x = (double*)malloc(n*sizeof(double));
  a = (double*)malloc(n*sizeof(double));

  srand(5489);

  int reps = 2;

  int rep=0;
  for(rep=0; rep<reps; ++rep) {
    int i=0;
    x[0] = 0;
    for(i=0; i<n;++i) {
      a[i] = i*i*i + (rand() % 10) * i;
      if(i)
        x[i] = x[i-1] + (rand() % 100) / 100.;
    }

    test_d(n, 3, x, a);
  }

  free(x);
  free(a);

  return 0;
}

void test_d(int n, int k, double *x, double *a)
{
  double *b, *b1;
  cs *D, *Dt;

  b = (double*)calloc(n, sizeof(double));
  b1 = (double*)calloc(n-k,sizeof(double));

  printf("x-------------------------------\n");
  print_array(x,n);

  printf("a-------------------------------\n");
  print_array(a,n);

  tf_dx(x,n,k,a,b);
  printf("D^k*a (in place)-----------------\n");
  print_array(b,n-k);

  printf("D^%d-------------------------\n", k);
  D = tf_calc_dk(n, k-1, x); /* D^k, not D^(k-1) */
  /* cs_print(D, 0); */
  cs_gaxpy(D,a,b1);
  printf("D^k*a (multiplication)-----------\n");
  print_array(b1,n-k);

  double err = max_diff(b, b1, n-k);
  printf("Dx err=%f\n", err);
  if(err > 1e-10)
    printf("D Test failed\n");


  tf_dtx(x,n,k,a,b);
  printf("(D^k)'*a (in place)-----------------\n");
  print_array(b,n);

  free(b1);
  b1 = (double*)calloc(n,sizeof(double));
  Dt = cs_transpose(D,1);
  cs_gaxpy(Dt,a,b1);
  printf("(D^k)'*a (multiplication)-----------\n");
  print_array(b1,n);

  err = max_diff(b, b1, n);
  printf("Dtx err=%f\n", err);
  if(err > 1e-10)
    printf("Dt Test failed\n");

  free(b);
  free(b1);
  cs_spfree(D);

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
