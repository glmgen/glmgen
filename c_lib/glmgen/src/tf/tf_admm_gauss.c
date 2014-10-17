#include "tf.h"

void tf_admm_gauss (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol,
       gqr * sparseQR)
{
  
  /* Some things that will be useful during our iterations */
  double *v = (double*)malloc(n*sizeof(double));
  double *z = (double*)malloc(n*sizeof(double));

  double pobj, loss, pen;
  int i, verb, it;

  verb = 0;
  if (verb) printf("Iteration\tObjective\tPenalty\n");

  for(it=0; it < max_iter; it++)
  {
    for (i=0; i < n-k; i++)
    {
      v[i] = alpha[i] + u[i];
    }
    /* z = \tilde{D}^T v = k D^T \Delta^-1 v */
    if( k > 0 )
      for(i=0; i < n-k; i++)
      {
        v[i] = v[i] * k/( x[k+i] - x[i] );
      }
    tf_dtx(x,n,k,v,z);

    for (i=0; i<n; i++)
    {
      /* beta[i] = w[i]*y[i] + rho*z[i]; */
      beta[i] = y[i] + rho*z[i]; /* Wy, i.e, y has weights */
    }
    /* Solve the least squares problem with sparse QR */
    glmgen_qrsol(sparseQR, beta);

    /* Update alpha: 1d fused lasso
     * Build the response vector */
    tf_dx(x,n,k,beta,v);
    if( k > 0 )
      for(i=0; i < n-k; i++)
      {
        v[i] = v[i] * k/( x[k+i] - x[i] );
      }
    for (i=0; i<n-k; i++)
    {
      z[i] = v[i]-u[i];
    }
    /* Use Nick's DP algorithm */
    tf_dp(n-k,z,lam/rho,alpha);

    /* Update u: dual update */
    for (i=0; i<n-k; i++)
    {
      u[i] = u[i]+alpha[i]-v[i];
    }

    /* Compute objective, if we are told to
     * Compute loss */
    loss = 0;
    for (i=0; i<n; i++)
    {
      /* loss += w[i]*(y[i]-beta[i])*(y[i]-beta[i]); */
      if( !( fabs(w[i]) < 1e-14 ) )
        loss += (y[i]-w[i]*beta[i])*(y[i]-w[i]*beta[i])/w[i];
    }
    /* Compute penalty */
    tf_dx(x,n,k+1,beta,z); /* IMPORTANT: use k+1 here! */
    pen = 0;
    for (i=0; i<n-k-1; i++)
    {
      pen += fabs(z[i]);
    }
    pobj = loss/2+lam*pen;
    obj[it] = pobj;

    if (verb) printf("%i\t%0.5e\t%0.5e\n",it,pobj, pen);

    /* Stop if relative difference of objective values <= obj_tol */
    if(it > 0)
    {
      if( fabs(pobj - obj[it-1]) < fabs(pobj) * obj_tol )
      {
        break;
      }
    }
  }

  *iter = it;

  free(v);
  free(z);
}


