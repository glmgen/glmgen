#include "tf.h"
#include <math.h>
#include <float.h>

/* Dynamic programming algorithm for the 1d fused lasso problem
   (Ryan's implementation of Nick Johnson's algorithm) */
void tf_dp (int n, double *y, double lam, double *beta)
{
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  int afirst;
  int alast;
  double bfirst;
  double blast;
  int alo;
  double blo;
  int ahi;
  double bhi;
  double * x;
  int * a;
  double * b;
  double *tm;
  double *tp;

  /* Take care of a few trivial cases */
  if (n==0) return;
  double maxy = fabs(y[0]);
  for(i=1; i < n; i++) maxy = max(maxy, fabs(y[i]));

  if (n==1 || lam <= DBL_EPSILON * 100*maxy)
  {
    for (i=0; i<n; i++) beta[i] = y[i];
    return;
  }
  

  x = (double*) malloc(2*n*sizeof(double));
  a = (int*) malloc(2*n*sizeof(int));
  b = (double*) malloc(2*n*sizeof(double));

  /* These are the knots of the back-pointers */
  tm = (double*) malloc((n-1)*sizeof(double));
  tp = (double*) malloc((n-1)*sizeof(double));

  /* We step through the first iteration manually */
  tm[0] = -lam+y[0];
  tp[0] = lam+y[0];
  l = n-1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -y[0]+lam;
  a[r] = -1;
  b[r] = y[0]+lam;
  afirst = 1;
  bfirst = -lam-y[1];
  alast = -1;
  blast = -lam+y[1];

  /* Now iterations 2 through n-1 */

  for (k=1; k<n-1; k++)
  {
    /* Compute lo: step up from l until the
       derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++)
    {
      if (alo*x[lo]+blo > -lam) break;
      alo += a[lo];
      blo += b[lo];
    }
   
    if(alo == 0) printf("1 k=%d\t  alo == 0\n",k);
    /* Compute the negative knot */
    tm[k] = (-lam-blo)/alo;
    l = lo-1;
    x[l] = tm[k];

    /* Compute hi: step down from r until the
       derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi=r; hi>=l; hi--)
    {
      if (-ahi*x[hi]-bhi < lam ) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    if(ahi == 0) printf("2 k=%d\t  ahi == 0, hi=%d, l=%d\n",k, hi,l);
    /* Compute the positive knot */
    tp[k] = (lam+bhi)/(-ahi);
    r = hi+1;
    x[r] = tp[k];

    if(alo == 0 || ahi == 0) printf("lam=%g, maxy=%g, lam/maxy=%g\n", lam, maxy, lam/maxy);
    /* Update a and b */
    a[l] = alo;
    b[l] = blo+lam;
    a[r] = ahi;
    b[r] = bhi+lam;
    afirst = 1;
    bfirst = -lam-y[k+1];
    alast = -1;
    blast = -lam+y[k+1];
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo=l; lo<=r; lo++)
  {
    if (alo*x[lo]+blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n-1] = -blo/alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (k=n-2; k>=0; k--)
  {
    if (beta[k+1]>tp[k]) beta[k] = tp[k];
    else if (beta[k+1]<tm[k]) beta[k] = tm[k];
    else beta[k] = beta[k+1];
  }

  /* Done! Free up memory */
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
}
