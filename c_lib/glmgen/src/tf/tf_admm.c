#include "tf.h"
#include "utils.h"

void tf_admm (double * y, double * x, double * w, int n, int k, int family,
              int max_iter, int lam_flag, double * lambda,
              int nlambda, double lambda_min_ratio, double * beta,
              double * obj, int * iter, int * status, double rho, double obj_tol)
{
  int i;
  double max_lam;
  double min_lam;
  double * temp_n;
  double * beta_max;
  double * alpha;
  double * u;
  double * resid;

  cs * D;
  cs * Dt;
  cs * Dk;
  cs * Dkt;
  cs * DktDk;
  gqr * Dt_qr;
  gqr * Dkt_qr;

  beta_max = (double *) malloc(n * sizeof(double));
  temp_n = (double *) malloc(n * sizeof(double));
  alpha = (double *) malloc(n * sizeof(double)); /* only size n-k, but we use extra buffer */
  u = (double *) malloc(n * sizeof(double)); /* only size n-k, but we use extra buffer */

  /* Assume w does not have zeros */
  for(i = 0; i < n; i++) temp_n[i] = 1/sqrt(w[i]);

  D = tf_calc_dk(n, k+1, x);
  Dk = tf_calc_dktil(n, k, x);
  Dt = cs_transpose(D, 1);
  diag_times_sparse(Dt, temp_n); /* Dt = W^{-1/2} Dt */
  Dkt = cs_transpose(Dk, 1);
  Dt_qr = glmgen_qr(Dt);
  Dkt_qr = glmgen_qr(Dkt);
  DktDk = cs_multiply(Dkt,Dk); 

  /* Determine the maximum lambda in the path, and initiate the path if needed
     using the input lambda_min_ratio and equally spaced log points. */

  max_lam = tf_maxlam(n, y, Dt_qr, w, family);
  if (!lam_flag)
  {
    min_lam = max_lam * lambda_min_ratio;
    lambda[0] = max_lam;
    for (i = 1; i < nlambda; i++)
      lambda[i] = exp((log(max_lam) * (nlambda - i) + log(min_lam) * i) / nlambda);
  }

  rho = rho * pow( (x[n-1] - x[0])/n, (double)k);

  /* Initiate alpha and u for a warm start */
  if (lambda[0] < max_lam * 1e-5)
  {
    for (i = 0; i < n - k; i++)
    {
      alpha[i] = 0; 
      u[i] = 0; 
    }
  } else {

    /* beta_max */
    for (i = 0; i < n; i++) temp_n[i] = -sqrt(w[i]) * y[i];
    glmgen_qrsol (Dt_qr, temp_n);    
    for (i = 0; i < n; i++) beta_max[i] = 0;
    cs_gaxpy(Dt, temp_n, beta_max); 
    // Dt has a W^{-1/2}, so in the next step divide by sqrt(w) instead of w.
    for (i = 0; i < n; i++) beta_max[i] = y[i] - beta_max[i]/sqrt(w[i]);

    /* alpha_max */
    tf_dxtil(x, n, k, beta_max, alpha);

    /* u_max */
    double exp_i;
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        for (i = 0; i < n; i++) u[i] = w[i] * (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;

      case FAMILY_LOGISTIC:
        for (i = 0; i < n; i++) {
          exp_i = exp(beta_max[i]);
          u[i] = exp_i / ((1+exp_i)*(1+exp_i)) * w[i] * (beta_max[i] - y[i]) / (rho * lambda[0]);
        }
        break;

      case FAMILY_POISSON:
        for (i = 0; i < n; i++) u[i] = exp(beta_max[i]) * w[i] *(beta_max[i] - y[i]) / (rho * lambda[0]);
        break;
    }
    glmgen_qrsol (Dkt_qr, u);
  }

  int j;

  /* Iterate lower level functions over all lambda values;
     the alpha and u vectors get used each time of subsequent
     warm starts */

  for (i = 0; i < nlambda; i++)
  {
    /* warm start */
    double * beta_init = (i == 0) ? beta_max : beta + (i-1)*n;
    for(j = 0; j < n; j++) beta[i*n + j] = beta_init[j];      
    
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        tf_admm_gauss(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
                      u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol, DktDk);
        break;

      case FAMILY_LOGISTIC:
        tf_admm_logistic(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
			 u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol, DktDk);
        break;

      case FAMILY_POISSON:
        tf_admm_poisson(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
		        u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol, DktDk);
        break;
    }
    
    /* If there any NaNs in beta: reset beta, alpha, u */
    if(!has_no_nan(beta + i * n, n)) {
      for(j = 0; j < n; j++) beta[i*n + j] = 0;
      for(j = 0; j < n-k; j++) { alpha[j] = 0; u[j] = 0; }
      status[i] = 1;
    }
  }

  cs_spfree(D);
  cs_spfree(Dt);
  cs_spfree(Dk);
  cs_spfree(Dkt);
  cs_spfree(DktDk);
  glmgen_gqr_free(Dt_qr);
  glmgen_gqr_free(Dkt_qr);

  free(temp_n);
  free(beta_max);
  free(alpha);
  free(u);
}
