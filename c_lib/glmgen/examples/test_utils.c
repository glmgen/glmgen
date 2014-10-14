#include "test_utils.h"
#include "tf.h"
#include "cs.h"
#include <math.h>

double rel_diff(double x, double y)
{
  double max_abs = fabs(x) > fabs(y) ? fabs(x) : fabs(y);

  return max_abs == 0 ? 0 : fabs(x-y)/max_abs;
}
double max_diff(double *x, double *y, int n)
{
  int i;
  double e=0;
  double c;
  for(i=0; i<n; i++)
  {
    c = rel_diff(x[i], y[i]);
    e = c > e ? c : e;
  }

  return e;
}

void print_array(double *x, int n)
{
  int i;

  for(i=0; i<n; i++)
  {
    printf("%f\n", x[i]);
  }
}

double norm(double * x, int n)
{
  if(!x)
    return 0;
  double sum = 0.;
  int i=0;
  for(i=0; i<n; ++i)
    sum += (x[i] * x[i]);

  return sqrt(sum);
}
double residual(cs * A, double * x, double * b)
{
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

