#include "tf.h"

void tf_admm_gauss (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, int vary_rho, double obj_tol,
       gqr * sparseQR, cs * DktDk)
{
  
  /* Some things that will be useful during our iterations */
  double *v = (double*)malloc(n*sizeof(double));
  double *z = (double*)malloc(n*sizeof(double));
  double *alpha_old = (double*)calloc(n, sizeof(double));
  gqr * qrfac;
  int new_qrfac;
  

  double pobj, loss, pen;
  int i, verb, it;

  double rho_mult;
  double rho_fac = 5.;
  double mu = 10.;

  cs * kernmat;
  double primal_res, dual_res;

  qrfac = sparseQR;
  new_qrfac = 0;

  verb = 0;
  if (verb) printf("Iteration\tObjective\tPenalty\n");

  for(it=0; it < max_iter; it++)
  {
    for (i=0; i < n-k; i++)
    {
      v[i] = alpha[i] + u[i];
    }

    tf_dtxtil(x,n,k,v,z);

    for (i=0; i<n; i++)
    {
      /* beta[i] = w[i]*y[i] + rho*z[i]; */
      beta[i] = y[i] + rho*z[i]; /* Wy, i.e, y has weights */
    }
    /* Solve the least squares problem with sparse QR */
    glmgen_qrsol(qrfac, beta);

    int num_nans;
    if(verb) {
      printf("it = %d\tb = %g", it, l2norm(beta,n));
      num_nans = count_nans(beta,n); if(num_nans) printf(", #nans=%d", num_nans);
    }

    /* Update alpha: 1d fused lasso
     * Build the response vector */
    tf_dxtil(x,n,k,beta,v);
    for (i=0; i<n-k; i++)
    {
      z[i] = v[i]-u[i];
    }
    memcpy(alpha_old, alpha, (n-k) * sizeof(double));
    /* Use Nick's DP algorithm */
    tf_dp(n-k,z,lam/rho,alpha);

    if(verb) {
      printf("\taold = %g\t z=%g", l2norm(alpha_old,n-k), l2norm(z,n-k));
      num_nans = count_nans(alpha,n); if(num_nans) printf(", #nans=%d", num_nans);
    }    
    /* Update u: dual update */
    for (i=0; i<n-k; i++)
    {
      u[i] = u[i]+alpha[i]-v[i];
    }
    if(verb) {
      printf("\tu = %g", l2norm(u,n-k));
      num_nans = count_nans(u,n); if(num_nans) printf(", #nans=%d", num_nans);
    }    
    if(verb) printf("\n");
    /* Compute objective, if we are told to
     * Compute loss */
    loss = 0;
    for (i=0; i<n; i++)
    {
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

    /*    if (verb) printf("%i\t%0.5e\t%0.5e\n",it,pobj, pen);*/

    /* Stop if relative difference of objective values <= obj_tol */
    if(it > 0)
    {
      if( fabs(pobj - obj[it-1]) < fabs(pobj) * obj_tol )
      {
        break;
      }
    }
    if( vary_rho ) {

      /*printf("v=%g\t alpha=%g, ", l2norm(v,n-k), l2norm(alpha,n-k));*/
      /* Compute primal residual \tilde{D}^{(k)} \beta - \alpha */
      for(i=0; i < n-k; i++) v[i] = v[i] - alpha[i];
      primal_res = l2norm(v,n-k);    

      /* Compute dual residual */
      for(i=0; i < n-k; i++) alpha_old[i] = alpha[i] - alpha_old[i];
      tf_dtxtil(x,n,k,alpha_old,alpha_old);
      dual_res = rho * l2norm(alpha_old, n);

      /* printf("r=%g\t s=%g\t rho=%g\n", primal_res, dual_res, rho); */

      /* Update rho, QR factorization and scale u */
      rho_mult = 1;
      if(primal_res > mu * dual_res)
        rho_mult = rho_fac;
      else if(dual_res > mu * primal_res)
        rho_mult = 1. / rho_fac;

      rho *= rho_mult;

      if(rho_mult != 1){
        if(new_qrfac) glmgen_gqr_free(qrfac);
        kernmat = scalar_plus_diag(DktDk, rho, w);
        qrfac = glmgen_qr(kernmat);
        new_qrfac = 1;
        cs_spfree(kernmat);

        /* scale u */
        for(i=0; i < n; i++) u[i] /= rho_mult;

      }
    }
  }

  *iter = it;

  free(v);
  free(z);
  free(alpha_old);
  if(new_qrfac) glmgen_gqr_free(qrfac);
}


