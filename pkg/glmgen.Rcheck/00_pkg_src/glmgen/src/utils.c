#include "utils.h"

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

  for (i=0; i<n; i++) {
    if (y[i]>lam) beta[i] = y[i]-lam;
    else if (y[i]<-lam) beta[i] = y[i]+lam;
    else beta[i]=0;
  }
}
