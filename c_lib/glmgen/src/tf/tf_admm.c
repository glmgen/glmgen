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
  double * u;
  double * Dy;
  double * resid;

  cs * D;
  cs * Dt;
  cs * DDt;
  cs * Dk;
  cs * Dkt;
  cs * DkDkt;
  cs * DktDk;
  gqr * DDt_qr;
  gqr * DkDkt_qr;
  gqr * DktDk_qr;

  cs * kernmat;
  gqr * kernmat_qr;

  beta_max = (double *) malloc(n * sizeof(double));
  temp_n = (double *) malloc(n * sizeof(double));
  alpha = (double *) malloc(n * sizeof(double)); /* only size n-k, but may need extra buffer */
  u = (double *) malloc((n - k) * sizeof(double));
  Dy = (double *) malloc((n - k - 1) * sizeof(double));
  resid = (double *) malloc(n * sizeof(double));

  D = tf_calc_dk(n, k+1, x);
  Dt = cs_transpose(D, 1);
  DDt = cs_multiply(D,Dt);
  Dk = tf_calc_dktil(n, k, x);
  Dkt = cs_transpose(Dk, 1);
  DkDkt = cs_multiply(Dk,Dkt);
  DktDk = cs_multiply(Dkt,Dk);
  DDt_qr = glmgen_qr(DDt);
  DkDkt_qr = glmgen_qr(DkDkt);
  DktDk_qr = glmgen_qr(DktDk);


  for (i = 0; i < n - k - 1; i++) Dy[i] = 0;
  cs_gaxpy(D, y, Dy);

  /* Determine the maximum lambda in the path, and initiate the path if needed
     using the input lambda_min_ratio and equally spaced log points. */
  max_lam = ts_maxlam(n - k-1, Dy, DDt_qr, family);
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
    for (i = 0; i < n-k-1; i++) temp_n[i] = -1*Dy[i];
    for (i = 0; i < n; i++) beta_max[i] = y[i];
    glmgen_qrsol (DDt_qr, temp_n);
    cs_gaxpy(Dt, temp_n, beta_max);

    /* alpha_max */
    tf_dx(x, n, k, beta_max, alpha);

    /* u_max */
    double exp_i;
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        for (i = 0; i < n; i++) temp_n[i] = (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;

      case FAMILY_LOGISTIC:
        for (i = 0; i < n; i++) {
          exp_i = exp(beta_max[i]);
          temp_n[i] = exp_i / ((1+exp_i)*(1+exp_i)) * (beta_max[i] - y[i]) / (rho * lambda[0]);
        }
        break;

      case FAMILY_POISSON:
        for (i = 0; i < n; i++) temp_n[i] = exp(beta_max[i]) * (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;
    }
    for (i = 0; i < n - k; i++) u[i] = 0;
    cs_gaxpy(Dk, temp_n, u);
    glmgen_qrsol (DkDkt_qr, u);

  }


  int j;
  /* Iterate lower level functions over all lambda values;
     the alpha and u vectors get used each time of subsequent
     warm starts */
  for (i = 0; i < nlambda; i++)
  {
    /* warm start */
    if(i == 0)
      for(j = 0; j < n; j++) beta[j] = beta_max[j];
    else
      for(j = 0; j < n; j++) beta[i*n + j] = beta[(i-1)*n + j];

    switch (family)
    {
      case FAMILY_GAUSSIAN:
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
  cs_spfree(DDt);
  cs_spfree(Dk);
  cs_spfree(Dkt);
  cs_spfree(DkDkt);
  cs_spfree(DktDk);
  glmgen_gqr_free(DDt_qr);
  glmgen_gqr_free(DkDkt_qr);
  glmgen_gqr_free(DktDk_qr);

  free(temp_n);
  free(beta_max);
  free(alpha);
  free(u);
  free(Dy);
  free(resid);
}
