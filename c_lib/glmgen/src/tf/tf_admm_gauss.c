#include "tf.h"

void tf_admm_gauss (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, cs * DktDk)
{  
  /* Special case for k=0: skip the ADMM algorithm */
  if (k==0) 
  {
    /* Use Nick's DP algorithm */
    tf_dp(n,y,lam,beta);

    int i;
    double loss, pen;
    double *db = (double*)malloc(n*sizeof(double));

    /* Compute objective */
    loss = 0; pen = 0;
    for (i=0; i<n; i++) loss += w[i]*(y[i]-beta[i])*(y[i]-beta[i]);
    tf_dx(x,n,k+1,beta,db); /* IMPORTANT: use k+1 here! */
    for (i=0; i<n-k-1; i++) pen += fabs(db[i]);
    obj[0] = loss/2+lam*pen;

    free(db);
    return;
  }

  /* Otherwise we run our ADMM routine */

  /* Construct the kernel matrix and its QR decomposition */
  cs * kernmat = scalar_plus_diag(DktDk, rho, w);
  gqr * kernmat_qr = glmgen_qr(kernmat);
  
  /* Other variables that will be useful during our iterations */
  double *v = (double*)malloc(n*sizeof(double));
  double *z = (double*)malloc(n*sizeof(double));
  double loss, pen;
  int i, verb, it;

  verb = 0;
  if (verb) printf("Iteration\tObjective\tPenalty\n");

  for(it=0; it < max_iter; it++)
  {
    /* Update beta: banded linear system (kernel matrix) */
    for (i=0; i < n-k; i++) v[i] = alpha[i] + u[i];
    tf_dtxtil(x,n,k,v,z);
    for (i=0; i<n; i++) beta[i] = w[i]*y[i] + rho*z[i];
    /* Solve the least squares problem with sparse QR */
    glmgen_qrsol(kernmat_qr, beta);

    int num_nans;
    if(verb) {
      printf("it = %d\tb = %g", it, beta[0]);
      num_nans = count_nans(beta,n); if(num_nans) printf(", #nans=%d", num_nans);
    }
    
    /* Update alpha: 1d fused lasso
     * Build the response vector */
    tf_dxtil(x,n,k,beta,v);
    for (i=0; i<n-k; i++)
    {
      z[i] = v[i]-u[i];
    }
    /* Use Nick's DP algorithm */
    tf_dp(n-k,z,lam/rho,alpha);

    if(verb) {
      printf("\ta = %g", alpha[0]);
      num_nans = count_nans(alpha,n); if(num_nans) printf(", #nans=%d", num_nans);
    }    
    /* Update u: dual update */
    for (i=0; i<n-k; i++)
    {
      u[i] = u[i]+alpha[i]-v[i];
    }
    if(verb) {
      printf("\tu = %g", u[0]);
      num_nans = count_nans(u,n); if(num_nans) printf(", #nans=%d", num_nans);
    }    
    if(verb) printf("\n");

    /* Compute objective */
    /* Compute loss */
    loss = 0;
    for (i=0; i<n; i++)
    {
      loss += w[i]*(y[i]-beta[i])*(y[i]-beta[i]);
    }
    /* Compute penalty */
    tf_dx(x,n,k+1,beta,z); /* IMPORTANT: use k+1 here! */
    pen = 0;
    for (i=0; i<n-k-1; i++)
    {
      pen += fabs(z[i]);
    }
    obj[it] = loss/2+lam*pen;

    //if (verb) printf("%i\t%0.5e\t%0.5e\n",it,obj[it], pen);

    /* Stop if relative difference of objective values <= obj_tol */
    if(it > 0)
    {
      if( fabs(obj[it] - obj[it-1]) < fabs(obj[it]) * obj_tol ) break;
    }
  }

  *iter = it;

  cs_spfree(kernmat);
  glmgen_gqr_free(kernmat_qr);
  free(v);
  free(z);
}


