#include "tf.h"

void tf_admm_gauss (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR)
{
  /* Some things that will be useful during our iterations */
  double *w = (double*)malloc(n*sizeof(double));
  double *z = (double*)malloc(n*sizeof(double));
  int *iter = (int*)malloc(sizeof(double));

  double pobj, loss, pen;

  int verb = 1; /* TODO */
  if (verb) printf("Iteration\tObjective");

  for (*iter=1; *iter<=max_iter; (*iter)++)
  {
    /* Update beta: banded least squares
     * Build the response vector
     * TODO: First add alpha and u and then pre-multiply by dt */
    tf_dtx(x,n,k,alpha,w);
    tf_dtx(x,n,k,u,z);
    int i=0;
    for (i=0; i<n; i++)
    {
      beta[i] = y[i] + rho*(w[i]+z[i]);
    }
    /* Solve the least squares problem with sparse QR */
    glmgen_qrsol(sparseQR, beta);

    /* Update alpha: 1d fused lasso
     * Build the response vector */
    tf_dx(x,n,k,beta,w);
    for (i=0; i<n-k; i++)
    {
      z[i] = w[i]-u[i];
    }
    /* Use Nick's DP algorithm */
    tf_dp(n-k,z,lam/rho,alpha);

    /* Update u: dual update */
    for (i=0; i<n-k; i++)
    {
      u[i] = u[i]+alpha[i]-w[i];
    }

    /* Compute objective, if we are told to
     * Compute loss */
    loss = 0;
    for (i=0; i<n; i++)
    {
      loss += (y[i]-beta[i])*(y[i]-beta[i]);
    }
    /* Compute penalty */
    tf_dx(x,n,k+1,beta,z); /* IMPORTANT: use k+1 here! */
    pen = 0;
    for (i=0; i<n-k-1; i++)
    {
      pen += fabs(z[i]);
    }
    pobj = loss/2+lam*pen;
    obj[(*iter)-1] = pobj;

    if (verb) printf("%i\t%0.5e\n",*iter,pobj);

    /* NEW: figure out when to stop
     * Based on a relative difference of objective values
     * being <= obj_tol */
    if(*iter > 1)
    {
      if( pobj < obj_tol || (obj[(*iter)-1] - pobj) / pobj < obj_tol )
      {
        break;
      }
    }
  }

  /* Clip the iteration counter at max_iter */
  if (*iter>max_iter) *iter = max_iter;

  /* TODO: Can we avoid allocating these arrays in each iteration? */
  free(w);
  free(z);
}


