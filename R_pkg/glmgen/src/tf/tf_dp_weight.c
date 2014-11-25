#include "tf.h"

/* Dynamic programming algorithm for the weighted 1d fused lasso problem 
   (Ryan's implementation of Nick Johnson's algorithm) */
void tf_dp_weight (int n, double *y, double *w, double lam, double *beta) 
{
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double alast;
  double bfirst;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;
  double *x;
  double *a;
  double *b;
  double *tm;
  double *tp;

  /* Take care of a few trivial cases */
  if (n==0) return;
  if (n==1 || lam==0)
  {
    for (i=0; i<n; i++) beta[i] = y[i];
    return;
  }
  
  /* Now deal with zero weights  */

  x = (double*) malloc(2*n*sizeof(double));
  a = (double*) malloc(2*n*sizeof(double));
  b = (double*) malloc(2*n*sizeof(double));

  /* These are the knots of the back-pointers */
  tm = (double*) malloc((n-1)*sizeof(double));
  tp = (double*) malloc((n-1)*sizeof(double));

  /* We step through the first iteration manually */
  tm[0] = -lam/w[0]+y[0];
  tp[0] = lam/w[0]+y[0];
  l = n-1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = w[0];
  b[l] = -w[0]*y[0]+lam;
  a[r] = -w[0];
  b[r] = w[0]*y[0]+lam;
  afirst = w[1];
  bfirst = -w[1]*y[1]-lam;
  alast = -w[1];
  blast = w[1]*y[1]-lam;

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
   
   /* Compute hi: step down from r until the
       derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi=r; hi>=lo; hi--)
    {
      if (-ahi*x[hi]-bhi < lam) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lam-blo)/alo;
    l = lo-1;
    x[l] = tm[k];
 
    /* Compute the positive knot */
    tp[k] = (lam+bhi)/(-ahi);
    r = hi+1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo+lam;
    a[r] = ahi;
    b[r] = bhi+lam;
    afirst = w[k+1];
    bfirst = -w[k+1]*y[k+1]-lam;
    alast = -w[k+1];
    blast = w[k+1]*y[k+1]-lam;

    /* double check=0; */
    /* check += afirst+alast; */
    /* for (i=l; i<=r; i++) check += a[i]; */
    /* if (check!=0) printf("k=%i, check=%f\n",k,check); */
    /* if (alo==0) printf("k=%d, alo=%f, lo=%d, r=%d, tm[k]=%e, x[lo]=%e\n",k,alo,lo,r,tm[k],x[lo]); */
    /* if (ahi==0) printf("k=%d, ahi=%f, hi=%d, l=%d, tp[k]=%e, x[hi]=%e\n",k,ahi,hi,l,tp[k],x[hi]); */
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
