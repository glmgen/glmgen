#include "utils.h"
#include <math.h>

double choose(int n, int k)
{
  return 1;
}

int max(int a, int b)
{
  return ((a) > (b) ? (a) : (b));
}

int min(int a, int b)
{
  return ((a) < (b) ? (a) : (b));
}

void soft_thresh(int n, double *y, double lam, double *beta)
{
  int i;

  for (i=0; i<n; i++)
  {
    if (y[i]>lam) beta[i] = y[i]-lam;
    else if (y[i]<-lam) beta[i] = y[i]+lam;
    else beta[i]=0;
  }
}

double glmgen_factorial(int n)
{
  int i=0;
  double x=1;

  for(i=2; i<=n; i++)
  {
    x *= i;
  }
  return x;
}

int is_nan(double x) {
  return (x != x);
}

int count_nans(double * x, int n) {
  int i;
  int num_nans = 0;

  for(i=0; i < n; i++) {
    num_nans += is_nan(x[i]);
  }
  return num_nans;
}

int has_nan(double * x, int n) {
  int i;
  for(i=0; i < n; i++) {
    if(is_nan(x[i]))
      return 1;
  }
  return 0;
}

double l1norm(double * x, int n) {
  int i;
  double s;
  s = 0;
  for(i=0; i < n; i++)
    s += fabs(x[i]);

  return s;
}
double l2norm(double * x, int n) {
  int i;
  double s;
  s = 0;
  for(i=0; i < n; i++)
    s += x[i] * x[i];

  return sqrt(s);
}
