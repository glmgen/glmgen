#include "tf.h"

void tf_admm (double * y, double * x, double * w, int n, int k, int family,
              int max_iter, int lam_flag, int obj_flag,  double * lambda,
              int nlambda, double lambda_min_ratio, double * beta,
              double * obj, double rho, double obj_tol)
{
  int i;
  double max_lam;
  double min_lam;
  double * temp_n_k;
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
  cs * DktDk;
  gqr * DDt_qr;
  gqr * DkDkt_qr;

  cs * kernmat;
  gqr * kernmat_qr;

  beta_max = (double *) malloc(n * sizeof(double));
  temp_n_k = (double *) malloc((n - k) * sizeof(double));
  alpha = (double *) malloc((n - k) * sizeof(double));
  u = (double *) malloc((n - k) * sizeof(double));
  Dy = (double *) malloc((n - k - 1) * sizeof(double));
  resid = (double *) malloc(n * sizeof(double));

  D = tf_calc_dk(n, k+1, x);
  Dt = cs_transpose(D, 1);
  DDt = cs_multiply(D,Dt);
  Dk = tf_calc_dk(n, k, x);
  Dkt = cs_transpose(Dk, 1);
  DktDk = cs_multiply(Dk,Dkt);
  DDt_qr = glmgen_qr(DDt);
  DkDkt_qr = glmgen_qr(DktDk);

  for (i = 0; i < n - k - 1; i++) Dy[i] = 0;
  cs_gaxpy(D, y, Dy);

  /* Determine the maximum lambda in the path, and initiate the path if needed
     using the input lambda_min_ratio and equally spaced log points. */
  max_lam = ts_maxlam(n - k, Dy, DDt_qr, family);
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

  /* Initiate beta_max */
  for (i = 0; i < n - k; i++) temp_n_k[i] = -1*Dy[i];
  for (i = 0; i < n; i++) beta_max[i] = y[i];
  glmgen_qrsol (DDt_qr, beta_max);
  cs_gaxpy(Dt, temp_n_k, beta_max);

  /* Initiate alpha and u for a warm start */
  if (lambda[0] < sqrt(max_lam))
  {
    for (i = 0; i < n - k; i++)
    {
      alpha[i] = 0;
      u[i] = 0;
    }
  } else {

    /* alpha_max */
    for (i = 0; i < n - k; i++) alpha[i] = 0;
    cs_gaxpy(Dk, beta_max, alpha);

    /* u_max */
    double exp_i;
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        for (i = 0; i < n - k; i++) temp_n_k[i] = (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;

      case FAMILY_LOGISTIC:
        for (i = 0; i < n - k; i++) {
          exp_i = exp(beta_max[i]);
          temp_n_k[i] = exp_i / ((1+exp_i)*(1+exp_i)) * (beta_max[i] - y[i]) / (rho * lambda[0]);
        }
        break;

      case FAMILY_POISSON:
        for (i = 0; i < n - k; i++) temp_n_k[i] = exp(beta_max[i]) * (beta_max[i] - y[i]) / (rho * lambda[0]);
        break;
    }
    for (i = 0; i < n - k; i++) u[i] = 0;
    cs_gaxpy(Dk, temp_n_k, u);
    glmgen_qrsol (DkDkt_qr, u);

  }

  /* Iterate lower level functions over all lambda values;
     the alpha and u vectors get used each time of subsequent
     warm starts */
  for (i = 0; i < nlambda; i++)
  {
    kernmat = scalar_plus_eye(DktDk, rho * lambda[i]);
    kernmat_qr = glmgen_qr(kernmat);
    switch (family)
    {
      case FAMILY_GAUSSIAN:
        tf_admm_gauss(y, x, w, n, k, max_iter, lambda[i], beta+i*n, alpha,
                      u, obj+i*max_iter, rho * lambda[i], obj_tol,
                      kernmat_qr);
        break;

      case FAMILY_LOGISTIC:
        break;

      case FAMILY_POISSON:
        break;
    }
    cs_spfree(kernmat);
    glmgen_gqr_free(kernmat_qr);
  }

  cs_spfree(D);
  cs_spfree(Dt);
  cs_spfree(DDt);
  cs_spfree(Dk);
  cs_spfree(Dkt);
  cs_spfree(DktDk);
  glmgen_gqr_free(DDt_qr);
  glmgen_gqr_free(DkDkt_qr);

  free(temp_n_k);
  free(beta_max);
  free(alpha);
  free(u);
  free(Dy);
  free(resid);
}