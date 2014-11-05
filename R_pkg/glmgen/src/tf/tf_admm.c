#include "tf.h"

#define TEST_FLAG 1

void tf_admm (double * y, double * x, double * w, int n, int k, int family,
              int max_iter, int lam_flag, int obj_flag,  double * lambda,
              int nlambda, double lambda_min_ratio, double * beta,
              double * obj, int * iter, double rho, double obj_tol)
{
  int i;
  double max_lam;
  double min_lam;
  double * temp_n;
  double * beta_max;
  double * alpha;
  double * alpha_max;
  double * u;
  double * u_max;
  double * resid;

  cs * D;
  cs * Dt;
  cs * Dk;
  cs * Dkt;
  cs * DktDk;
  gqr * Dt_qr;
  gqr * Dkt_qr;

  cs * kernmat;
  gqr * kernmat_qr;

  beta_max = (double *) malloc(n * sizeof(double));
  temp_n = (double *) malloc(n * sizeof(double));
  alpha = (double *) malloc(n * sizeof(double)); /* only size n-k, but may need extra buffer */
  u = (double *) malloc(n * sizeof(double));/* only size n-k, but may need extra buffer */
  alpha_max = (double *) malloc(n * sizeof(double));
  u_max = (double *) malloc(n * sizeof(double));

  resid = (double *) malloc(n * sizeof(double));

  D = tf_calc_dk(n, k+1, x);
  Dt = cs_transpose(D, 1);
  Dk = tf_calc_dktil(n, k, x);
  Dkt = cs_transpose(Dk, 1);
  DktDk = cs_multiply(Dkt,Dk);
  Dt_qr = glmgen_qr(Dt);
  Dkt_qr = glmgen_qr(Dkt);

  /* Determine the maximum lambda in the path, and initiate the path if needed
     using the input lambda_min_ratio and equally spaced log points. */
  max_lam = ts_maxlam(n, y, Dt_qr, family);
  if (!lam_flag)
  {
    min_lam = max_lam * lambda_min_ratio;
    lambda[0] = max_lam;
    if (nlambda != 1)
    {
      for (i = 1; i < nlambda; i++)
        lambda[i] = exp((log(max_lam) * (nlambda - i) + log(min_lam) * i) / nlambda);
    }
  }
  
  rho = rho * pow( (x[n-1] - x[0])/n, (double)k);

  /* Initiate alpha and u for a warm start */
  if (lambda[0] < max_lam * 1e-5)
  {
    for (i = 0; i < n - k; i++)
    {
      alpha[i] = 0; alpha_max[i] = 0;
      u[i] = 0; u_max[i] =0;
    }
  } else {

    /* beta_max */
    for (i = 0; i < n; i++) temp_n[i] = -y[i];
    for (i = 0; i < n; i++) beta_max[i] = y[i];
    glmgen_qrsol (Dt_qr, temp_n);
    cs_gaxpy(Dt, temp_n, beta_max);

    /* alpha_max */
    tf_dxtil(x, n, k, beta_max, alpha_max);
    memcpy(alpha, alpha_max, n*sizeof(double));

    /* u_max */
    double exp_i;
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        for (i = 0; i < n; i++) u_max[i] = (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;

      case FAMILY_LOGISTIC:
        for (i = 0; i < n; i++) {
          exp_i = exp(beta_max[i]);
          u_max[i] = exp_i / ((1+exp_i)*(1+exp_i)) * (beta_max[i] - y[i]) / (rho * lambda[0]);
        }
        break;

      case FAMILY_POISSON:
        for (i = 0; i < n; i++) u_max[i] = exp(beta_max[i]) * (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;
    }
    glmgen_qrsol (Dkt_qr, u_max);

    memcpy(u, u_max, n*sizeof(double));
  }

  int warm_start;
  int j;
  /* Iterate lower level functions over all lambda values;
     the alpha and u vectors get used each time of subsequent
     warm starts */
  for (i = 0; i < nlambda; i++)
  {
    /* warm start or cold start */
    double * beta_init = (i == 0) ? beta_max : beta + (i-1)*n;    
    
    warm_start = has_no_nan(beta_init, n);
    if(warm_start) {
      for(j = 0; j < n; j++) beta[i*n + j] = beta_init[j];      
    }
    else {
      for(j = 0; j < n; j++) beta[i*n + j] = beta_max[j];
      memcpy(alpha, alpha_max, n*sizeof(double));
      memcpy(u, u_max, n*sizeof(double));
    }
    
/*    if(i == 0)*/
/*      for(j = 0; j < n; j++) beta[j] = beta_max[j];*/
/*    else*/
/*      for(j = 0; j < n; j++) beta[i*n + j] = beta[(i-1)*n + j];*/
    
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        if(k == 0) {
          tf_dp(n, y, rho * lambda[i], beta+i*n);
          break;
        }
        kernmat = scalar_plus_eye(DktDk, rho * lambda[i]);
        kernmat_qr = glmgen_qr(kernmat);
        
        tf_admm_gauss(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
                      u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol,
                      kernmat_qr);
                      
        cs_spfree(kernmat);
        glmgen_gqr_free(kernmat_qr);
        break;

      case FAMILY_LOGISTIC:
        tf_admm_logistic(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
                      u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol);
        break;

      case FAMILY_POISSON:
        tf_admm_pois(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
                      u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol);
        break;
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
  free(alpha_max);
  free(u_max); 
  free(resid);
}
