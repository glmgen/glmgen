#include <math.h>
#include <stdlib.h>

#include "tf.h"
#include "matcomp.h"
#include "utils.h"
#include "int_codes.h"

void tf_admm (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj,
              double rho, double obj_tol)
{
  int i;
  int j;
  double max_lam;
  double min_lam;
  double * beta_max;
  double * dtd;
  double * u;
  double * alpha;
  csn * sparseQR;

  beta_max = (double *) malloc(n * sizeof(double));
  dtd = (double *) malloc(n * sizeof(double)); // TODO: what is the size here?
  u = (double *) malloc(n * sizeof(double));
  alpha = (double *) malloc(n * sizeof(double));

  max_lam = ts_maxlam(y, x, n, k, beta_max);
  if (!lam_flag)
  {
    min_lam = max_lam * lambda_min_ratio;
    for (i = 0; i < nlambda; i++)
      lambda[i] = exp((log(max_lam) * (nlambda - i) + log(min_lam) * i) / nlambda);
    lambda[0] = max_lam; // force this, for rounding errors
  }

  tf_calc_dtd(x, n, k, dtd);

  for(i = 0; i < nlambda; i++)
  {
    tf_getrho(&rho, lambda[i]);
    tf_calc_sparse_qr(n, k, rho, dtd, sparseQR);
    switch(family)
    {
      case FAMILY_GAUSSIAN:
        tf_admm_gauss(y, x, n, k, max_iter, lambda[i], beta+i*n, alpha, u, obj+i*max_iter, rho, obj_tol, sparseQR);
        break;

      case FAMILY_LOGISTIC:
        break;

      case FAMILY_POISSON:
        break;
    }
  }

  free(beta_max);
  free(dtd);
  free(u);
  free(alpha);
}

void tf_primal_dual (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj)
{

}

void tf_admm_gauss(double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR)
{

  // NEW: here we will have already computed a sparse QR
  // of the kernel matrix (rho*D^T*D + I). This will have been
  // done by our parent function
  // The actual sparseQR type will not be a (double *)

  // OLD: with banded Cholesky
  /* double *drow = (double*)malloc((k+1)*sizeof(double)); */
  /* builddrow(k,drow); */
  /* char uplo = 'L'; */
  /* int kd=k, nrhs=1, ldab=k+1, info=0; */
  /* double *dtd = (double*)malloc(ldab*n*sizeof(double)); */
  /* double *ab = (double*)malloc(ldab*n*sizeof(double)); */
  /* // Note: dtd and ab are stored in the packed representation for */
  /* // symmetric banded matrices, in *column-major order* (for use */
  /* // by LAPACK in Fortran!) */
  /* builddtd(n,k,drow,dtd); */
  /* buildab(n,k,dtd,ab,rho); */
  /* // Compute Cholesky factorization */
  /* F77_CALL(dpbtrf)(&uplo,&n,&kd,ab,&ldab,&info); */

  int i;
  int iter;
  double * w;
  double * z;
  double pobj;
  double loss;
  double pen;

  w = (double*)malloc(n*sizeof(double));
  z = (double*)malloc(n*sizeof(double));

  // printf("Iteration\tObjective");

  for (iter=1; iter<=max_iter; (iter)++)
  {
    // Update beta: banded least squares
    // Build the response vector
    dt(n, k, w, alpha, x);
    dt(n, k, z, u, x);
    for (i=0; i<n; i++)
    {
      beta[i] = y[i] + rho*(w[i]+z[i]);
    }
    // NEW: Solve the least squares problem with sparse QR
    // OLD: with banded Cholesky
    // F77_CALL(dpbtrs)(&uplo,&n,&k,&nrhs,ab,&ldab,beta,&n,&info);

    // Update alpha: 1d fused lasso
    // Build the response vector
    d(n, k, w, beta, x);
    for (i=0; i<n-k; i++)
    {
      z[i] = w[i]-u[i];
    }
    // Use Nick's DP algorithm
    tf_dp(n-k,z,lam/rho,alpha);

    // Update u: dual update
    for (i=0; i<n-k; i++)
    {
      u[i] = u[i]+alpha[i]-w[i];
    }

    // Compute objective, if we are told to
    // Compute loss
    loss = 0;
    for (i=0; i<n; i++)
    {
      loss += (y[i]-beta[i])*(y[i]-beta[i]);
    }
    // Compute penalty
    d(n, k+1, z, beta, x);// IMPORTANT: use k+1 here!
    pen = 0;
    for (i=0; i<n-k-1; i++)
    {
      pen += fabs(z[i]);
    }
    pobj = loss/2+lam*pen;
    obj[(iter)-1] = pobj;

    // printf("%i\t%0.5e\n",*iter,pobj);

    // NEW: figure out when to stop
    // Based on a relative difference of objective values
    // being <= obj_tol
  }

  // Clip the iteration counter at maxiter
  if (iter>max_iter) iter = max_iter;

  // NEW: free sparse QR
  // OLD: free Cholesky variables
  // free(dtd);
  // free(ab);

  free(w);
  free(z);
}

void tf_admm_logistic(double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR)
{

}

void tf_admm_pois(double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR)
{

}

// Dynamic programming algorithm for the 1d fused lasso problem
// (Ryan's implementation of Nick Johnson's algorithm)
void tf_dp(int n, double *y, double lam, double *beta)
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
  double * x;
  double * a;
  double * b;
  double *tm;
  double *tp;

  // Take care of a few trivial cases
  if (n==0) return;
  if (n==1 || lam==0)
  {
    for (i=0; i<n; i++) beta[i] = y[i];
    return;
  }

  x = (double*) malloc(2*n*sizeof(double));
  a = (double*) malloc(2*n*sizeof(double));
  b = (double*) malloc(2*n*sizeof(double));

  // These are the knots of the back-pointers
  tm = (double*) malloc((n-1)*sizeof(double));
  tp = (double*) malloc((n-1)*sizeof(double));

  // We step through the first iteration manually
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

  // Now iterations 2 through n-1

  for (k=1; k<n-1; k++)
  {
    // Compute lo: step up from l until the
    // derivative is greater than -lam
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++)
    {
      if (alo*x[lo]+blo > -lam) break;
      alo += a[lo];
      blo += b[lo];
    }

    // Compute the negative knot
    tm[k] = (-lam-blo)/alo;
    l = lo-1;
    x[l] = tm[k];

    // Compute hi: step down from r until the
    // derivative is less than lam
    ahi = alast;
    bhi = blast;
    for (hi=r; hi>=l; hi--)
    {
      if (-ahi*x[hi]-bhi < lam) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    // Compute the positive knot
    tp[k] = (lam+bhi)/(-ahi);
    r = hi+1;
    x[r] = tp[k];

    // Update a and b
    a[l] = alo;
    b[l] = blo+lam;
    a[r] = ahi;
    b[r] = bhi+lam;
    afirst = 1;
    bfirst = -lam-y[k+1];
    alast = -1;
    blast = -lam+y[k+1];
  }

  // Compute the last coefficient: this is where
  // the function has zero derivative
  alo = afirst;
  blo = bfirst;
  for (lo=l; lo<=r; lo++)
  {
    if (alo*x[lo]+blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n-1] = -blo/alo;

  // Compute the rest of the coefficients, by the
  // back-pointers
  for (k=n-2; k>=0; k--)
  {
    if (beta[k+1]>tp[k]) beta[k] = tp[k];
    else if (beta[k+1]<tm[k]) beta[k] = tm[k];
    else beta[k] = beta[k+1];
  }

  // Done! Free up memory
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
}

double ts_maxlam(double * y, double * x, int n, int k, double * beta_max)
{
  return 1;
}

void tf_calc_dtd(double * x, int n, int k, double * dtd)
{

}

void tf_getrho(double * rho, double lambda)
{
  rho[0] = lambda;
}

void tf_calc_sparse_qr(int n, int k, double rho, double * dtd, csn * sparseQR)
{

}

void tf_d(double *x, int n, int k,double *a, double *b)
{

}

void tf_dt(double *x, int n, int k, double *a, double *b)
{

}
