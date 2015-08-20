/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Veeranjaneyulu Sadhanala,           *
 *                       Ryan Tibshirani                                    *
 *                                                                          *
 * This file is part of the glmgen library / package.                       *
 *                                                                          *
 *   glmgen is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 2 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   glmgen is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with glmgen. If not, see <http://www.gnu.org/licenses/>. *
 ****************************************************************************/

/**
 * @file tf_admm.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-23
 * @brief Fitting trend filtering estimates with the ADMM algorithm.
 * Contains all of the functions specific to the ADMM algorithm
 * implementation; allows for Gaussian, binomial, and poisson
 * losses, as well as arbitrary sample weights.
 */

#include "tf.h"
#include "utils.h"

/**
 * @brief Default call to tf_admm.
 * Example of how to call tf_admm, taking only the response vector @p and
 * observation size @n as inputs. Users will likely need to adjust tf_admm_default
 * for their own needs, using it as a starting template.
 * @param y                    a vector of responses
 * @param n                    the length of y, i.e., the number of observations
 * @return  Returns a pointer to beta values, a column oriented (nlambda x n) array.
 * @note The user is responsible for freeing the array pointed to by the response.
 */
double * tf_admm_default(double * y, int n)
{
  int i;
  int j;

  double * x;
  double * w;
  int k;
  int family;
  int max_iter;
  int lam_flag;
  double * lambda;
  int nlambda;
  double lambda_min_ratio;
  int * df;
  double * beta;
  double * obj;
  int * iter;
  int * status;
  double rho;
  double obj_tol;
  double alpha_ls;
  double gamma_ls;
  int max_iter_ls;
  int max_iter_newton;
  int verbose;

  /* Set default constants */
  k = 2;
  family = FAMILY_GAUSSIAN;
  max_iter = 200;
  lam_flag = 0;
  nlambda = 50;
  lambda_min_ratio = 1e-5;
  rho = 1;
  obj_tol = 1e-4;
  alpha_ls = 0.5;
  gamma_ls = 0.8;
  max_iter_ls = 20;
  max_iter_newton = 30;
  verbose = 0;

  /* Allocate space for input arrays */
  x      = (double *) malloc(n * sizeof(double));
  w      = (double *) malloc(n * sizeof(double));
  lambda = (double *) malloc(nlambda * sizeof(double));
  df     = (int *)    malloc(nlambda * sizeof(int));
  beta   = (double *) malloc(nlambda * n * sizeof(double));
  obj    = (double *) malloc(nlambda * max_iter * sizeof(double));
  iter   = (int *)    malloc(nlambda * sizeof(int));
  status = (int *)    malloc(nlambda * sizeof(int));

  /* Fill x and w with default values */
  for (i = 0; i < n; i++)
  {
    x[i] = i;
    w[i] = 1;
  }

  /* Initalize output arrays with 0's */
  for (i = 0; i < nlambda; i++)
  {
    lambda[i] = 0;
    df[i] = 0;
    for (j = 0; j < n; j++) beta[i + j*nlambda] = 0;
    for (j = 0; j < max_iter; j++) obj[i + j*nlambda] = 0;
    iter[i] = 0;
    status[i] = 0;
  }

  tf_admm(x, y, w, n, k, family, max_iter, lam_flag, lambda,
          nlambda, lambda_min_ratio, df, beta, obj, iter,
          status, rho, obj_tol, alpha_ls, gamma_ls, max_iter_ls,
          max_iter_newton, verbose);

  /* Free allocated arrays (except beta; which is returned) */
  free(x);
  free(w);
  free(lambda);
  free(df);
  free(obj);
  free(iter);
  free(status);

  return(beta);
}

/**
 * @brief Main wrapper for fitting a trendfilter model.
 * Takes as input either a sequence of lambda tuning parameters, or the number
 * of desired lambda values. In the latter case the function will also calculate
 * a lambda sequence. The user must supply allocated memory to store the output,
 * with the function itself returning only @c void. For default values, and an
 * example of how to call the function, see the function tf_admm_default.
 *
 * @param x                    a vector of data locations; must be in increasing order
 * @param y                    a vector of responses
 * @param w                    a vector of sample weights
 * @param n                    the length of x, y, and w
 * @param k                    polynomial degree of the fitted trend; i.e., k=1 for linear
 * @param family               family code for the type of fit; family=0 for OLS
 * @param max_iter             maximum number of ADMM interations; ignored for k=0
 * @param lam_flag             0/1 flag for whether lambda sequence needs to be estimated
 * @param lambda               either a sequence of lambda when lam_flag=0, or empty
 *                             allocated space if lam_flag=1
 * @param nlambda              number of lambda values; need for both lam_flag=0 and 1
 * @param lambda_min_ratio     minimum ratio between min and max lambda; ignored for lam_flag=0
 * @param df                   allocated space of nlambda to store the output df values
 * @param beta                 allocated space of size n*nlambda to store the output coefficents
 * @param obj                  allocated space of size max_iter*nlambda to store the objective
 * @param iter                 allocated space of size nlambda to store the number of iterations
 * @param status               allocated space of size nlambda to store the status of each run
 * @param rho                  tuning parameter for the ADMM algorithm
 * @param obj_tol              stopping criteria tolerance
 * @param alpha_ls             for family != 0, line search tuning parameter
 * @param gamma_ls             for family != 0, line search tuning parameter
 * @param max_iter_ls          for family != 0, max number of iterations in line search
 * @param max_iter_newton      for family != 0, max number of iterations in inner ADMM
 * @param verbose              0/1 flag for printing progress
 * @return void
 * @see tf_admm_default
 */
void tf_admm (double * x, double * y, double * w, int n, int k, int family,
              int max_iter, int lam_flag, double * lambda,
              int nlambda, double lambda_min_ratio, int * df,
	      double * beta, double * obj, int * iter, int * status,
	      double rho, double obj_tol, double alpha_ls, double gamma_ls,
              int max_iter_ls, int max_iter_newton, int verbose)
{
  int i;
  int j;
  double max_lam;
  double min_lam;
  double * temp_n;
  double * beta_max;
  double * alpha;
  double * u;

  cs * D;
  cs * Dt;
  cs * Dk;
  cs * Dkt;
  cs * DktDk;
  gqr * Dt_qr;
  gqr * Dkt_qr;

  beta_max = (double *) malloc(n * sizeof(double));
  temp_n   = (double *) malloc(n * sizeof(double));
  alpha    = (double *) malloc(n * sizeof(double)); /* we use extra buffer (n vs n-k) */
  u        = (double *) malloc(n * sizeof(double)); /* we use extra buffer (n vs n-k) */

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
   * using the input lambda_min_ratio and equally spaced log points.
   */
  max_lam = tf_maxlam(n, y, Dt_qr, w);
  if (family == FAMILY_GAUSSIAN) {
    if (max_lam <= 1e-10 * l1norm(y,n)/n) {
      for (i=0; i<nlambda; i++) {
	for (j=0; j<n; j++) beta[i*n+j] = y[j];
	obj[i*max_iter] = 0;
	lambda[i] = 0;
      }
      // TODO free stuff
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
      return;
    }
  }
  else {
    max_lam += 1;
  }
  if (!lam_flag) {
    min_lam = max_lam * lambda_min_ratio;
    lambda[0] = max_lam;
    for (i = 1; i < nlambda; i++)
      lambda[i] = exp((log(max_lam) * (nlambda - i -1) + log(min_lam) * i) / (nlambda-1));
  }
  rho = rho * pow((x[n-1] - x[0])/n, (double)k);

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
    for (i = 0; i < n; i++) temp_n[i] = sqrt(w[i]) * y[i];
    glmgen_qrsol (Dt_qr, temp_n);
    for (i = 0; i < n; i++) beta_max[i] = 0;
    cs_gaxpy(Dt, temp_n, beta_max);
    /* Dt has a W^{-1/2}, so in the next step divide by sqrt(w) instead of w. */
    for (i = 0; i < n; i++) beta_max[i] = y[i] - beta_max[i]/sqrt(w[i]);

    /* Check if beta = weighted mean(y) is better than beta */
    double yc = 0, sumw = 0;
    for (i = 0; i < n; i++) yc += w[i] * y[i];
    for (i = 0; i < n; i++) sumw += w[i];
    yc /= sumw;
    switch(family) {
    case FAMILY_POISSON:
      yc = (yc > 0)? log(yc) : -DBL_MAX;
      break;
    case FAMILY_LOGISTIC:
      yc = (yc > 0) ? ( yc < 1 ? log(yc/(1-yc)) : DBL_MAX) : -DBL_MAX;
      break;
    default: break;
    }
    for (i = 0; i < n; i++) temp_n[i] = yc;
    double obj1 = tf_obj(x,y,w,n,k,max_lam,family,beta_max,alpha);
    double obj2 = tf_obj(x,y,w,n,k,max_lam,family,temp_n,alpha);
    if(obj2 < obj1) {
      for (i = 0; i < n; i++) beta_max[i] = yc;
    }

    /* alpha_max */
    tf_dxtil(x, n, k, beta_max, alpha);

    /* u_max */
    switch (family)
    {
    case FAMILY_GAUSSIAN:
      for (i = 0; i < n; i++) u[i] = w[i] * (beta_max[i] - y[i]) / (rho * lambda[0]);
      break;

    case FAMILY_LOGISTIC:
      for (i = 0; i < n; i++) {
        u[i] = logi_b2(beta_max[i]) * w[i] * (beta_max[i] - y[i]) / (rho * lambda[0]);
      }
      break;

    case FAMILY_POISSON:
      for (i = 0; i < n; i++) {
        u[i] = pois_b2(beta_max[i]) * w[i] *(beta_max[i] - y[i]) / (rho * lambda[0]);
      }
      break;

    default:
      for (i = 0; i < nlambda; i++) status[i] = 2;
      return;
    }

    glmgen_qrsol (Dkt_qr, u);
  }

  /* Iterate lower level functions over all lambda values;
   * the alpha and u vectors get used each time of subsequent
   * warm starts
   */
  for (i = 0; i < nlambda; i++)
  {
    /* warm start */
    double * beta_init = (i == 0) ? beta_max : beta + (i-1)*n;
    for(j = 0; j < n; j++) beta[i*n + j] = beta_init[j];

    switch (family)
    {
      case FAMILY_GAUSSIAN:
        tf_admm_gauss(x, y, w, n, k, max_iter, lambda[i], df+i, beta+i*n,
		      alpha, u, obj+i*max_iter, iter+i, rho * lambda[i],
		      obj_tol, DktDk, verbose);
        break;

      case FAMILY_LOGISTIC:
        tf_admm_glm(x, y, w, n, k, max_iter, lambda[i], df+i, beta+i*n,
		    alpha, u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol,
		    alpha_ls, gamma_ls, max_iter_ls, max_iter_newton,
                    DktDk, &logi_b, &logi_b1, &logi_b2, verbose);
        break;

      case FAMILY_POISSON:
        tf_admm_glm(x, y, w, n, k, max_iter, lambda[i], df+i, beta+i*n,
		    alpha, u, obj+i*max_iter, iter+i, rho * lambda[i], obj_tol,
		    alpha_ls, gamma_ls, max_iter_ls, max_iter_newton,
		    DktDk, &pois_b, &pois_b1, &pois_b2, verbose);
        break;
    }

    /* If there any NaNs in beta: reset beta, alpha, u */
    if (has_nan(beta + i*n, n))
    {
      double yc = 0, sumw = 0;
      for (i = 0; i < n; i++) yc += w[i] * y[i];
      for (i = 0; i < n; i++) sumw += w[i];
      yc /= sumw;
      switch(family) {
      case FAMILY_POISSON:
	yc = (yc > 0) ? log(yc) : -DBL_MAX;
	break;
      case FAMILY_LOGISTIC:
	yc = (yc > 0) ? ( yc < 1 ? log(yc/(1-yc)) : DBL_MAX) : -DBL_MAX;
	break;
      default: break;
      }
      for (j = 0; j < n; j++) beta[i*n + j] = yc;
      for (j = 0; j < n-k; j++) alpha[j] = 0;
      for (j = 0; j < n; j++) u[j] = w[j] * (beta[i*n+j] - y[j]) / (rho * lambda[i]);
      glmgen_qrsol (Dkt_qr, u);
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

/**
 * @brief Low level fitting routine for a Gaussian trend filtering problem.
 * Function used by tf_admm to fit a Gaussian ADMM trendfilter, or as a
 * subproblem by tf_admm_glm when using logistic or poisson losses. Fits
 * the solution for a single value of lambda. Most users will want to call
 * tf_admm, rather than tf_admm_gauss directly.
 *
 * @param x                    a vector of data locations; must be in increasing order
 * @param y                    a vector of responses
 * @param w                    a vector of sample weights
 * @param n                    the length of x, y, and w
 * @param k                    polynomial degree of the fitted trend; i.e., k=1 for linear
 * @param max_iter             maximum number of ADMM interations; ignored for k=0
 * @param lam                  the value of lambda
 * @param df                   allocated space for df value at the solution
 * @param beta                 allocated space for output coefficents; must pre-fill as it is used in warm start
 * @param alpha                allocated space for ADMM alpha variable; must pre-fill as it is used in warm start
 * @param u                    allocated space for ADMM u variable; must pre-fill as it is used in warm start
 * @param obj                  allocated space to store the objective; will fill at most max_iter elements
 * @param iter                 allocated space to store the number of iterations; will fill just one element
 * @param rho                  tuning parameter for the ADMM algorithm; set to 1 for default
 * @param obj_tol              stopping criteria tolerance; set to 1e-10 for default
 * @param DktDk                pointer to the inner product of DktDk
 * @param verbose              0/1 flag for printing progress
 * @return void
 * @see tf_admm
 */
void tf_admm_gauss (double * x, double * y, double * w, int n, int k,
		    int max_iter, double lam, int * df,
		    double * beta, double * alpha, double * u,
		    double * obj, int * iter,
		    double rho, double obj_tol, cs * DktDk, int verbose)
{
  int i;
  int d;
  int it;
  double *v;
  double *z;
  double *db;

  cs * kernmat;
  gqr * kernmat_qr;

  /* Special case for k=0: skip the ADMM algorithm */
  if (k==0)
  {
    /* Use Nick's DP algorithm, weighted version */
    tf_dp_weight(n,y,w,lam,beta);

    /* Compute df value */
    d = 1;
    for (i=0; i<n-1; i++) if (beta[i] != beta[i+1]) d += 1;
    *df = d;

    /* Compute objective */
    db = (double *) malloc(n*sizeof(double));
    obj[0] = tf_obj(x,y,w,n,k,lam,FAMILY_GAUSSIAN,beta,db);
    free(db);
    return;
  }

  /* Otherwise we run our ADMM routine */

  /* Construct the kernel matrix and its QR decomposition */
  kernmat = scalar_plus_diag(DktDk, rho, w);
  kernmat_qr = glmgen_qr(kernmat);

  /* Other variables that will be useful during our iterations */
  v = (double*) malloc(n*sizeof(double));
  z = (double*) malloc(n*sizeof(double));

  if (verbose) printf("\nlambda=%0.3e\n",lam);
  if (verbose) printf("Iteration\tObjective\n");

  for (it=0; it < max_iter; it++)
  {
    /* Update beta: banded linear system (kernel matrix) */
    for (i=0; i < n-k; i++) v[i] = alpha[i] + u[i];
    tf_dtxtil(x,n,k,v,z);
    for (i=0; i<n; i++) beta[i] = w[i]*y[i] + rho*z[i];
    /* Solve the least squares problem with sparse QR */
    glmgen_qrsol(kernmat_qr, beta);

    /* Update alpha: 1d fused lasso
     * Build the response vector */
    tf_dxtil(x,n,k,beta,v);
    for (i=0; i<n-k; i++) z[i] = v[i]-u[i];

    /* Use Nick's DP algorithm */
    tf_dp(n-k,z,lam/rho,alpha);

    /* Update u: dual update */
    for (i=0; i<n-k; i++) u[i] = u[i]+alpha[i]-v[i];

    /* Compute objective */
    obj[it] = tf_obj(x,y,w,n,k,lam,FAMILY_GAUSSIAN,beta,z);
    if (verbose) printf("%i\t%0.3e\n",it+1,obj[it]);

    /* Stop if relative difference of objective values < obj_tol */
    if (it > 0 && (fabs(obj[it] - obj[it-1]) < fabs(obj[it-1]) * obj_tol)) break;
  }

  *iter = it;

  /* Compute final df value, based on alpha */
  d = k+1;
  for (i=0; i<n-k; i++) if (alpha[i] != alpha[i+1]) d += 1;
  *df = d;

  cs_spfree(kernmat);
  glmgen_gqr_free(kernmat_qr);
  free(v);
  free(z);
}

/**
 * @brief Low level fitting routine for non-Gaussian trend filtering problems.
 * Can be configured to handle arbirary losses, as it takes the link function
 * and it first two derivaties as inputs. Fits the solution for a single value
 * of lambda. Most users will want to call tf_admm, rather than tf_admm_glm directly.
 *
 * @param x                    a vector of data locations; must be in increasing order
 * @param y                    a vector of responses
 * @param w                    a vector of sample weights
 * @param n                    the length of x, y, and w
 * @param k                    polynomial degree of the fitted trend; i.e., k=1 for linear
 * @param max_iter             maximum number of ADMM interations; ignored for k=0
 * @param lam                  the value of lambda
 * @param df                   allocated space for df value at the solution
 * @param beta                 allocated space for output coefficents; must pre-fill as it is used in warm start
 * @param alpha                allocated space for ADMM alpha variable; must pre-fill as it is used in warm start
 * @param u                    allocated space for ADMM u variable; must pre-fill as it is used in warm start
 * @param obj                  allocated space to store the objective; will fill at most max_iter elements
 * @param iter                 allocated space to store the number of iterations; will fill just one element
 * @param status               allocated space of size nlambda to store the status of each run
 * @param rho                  tuning parameter for the ADMM algorithm; set to 1 for default
 * @param obj_tol              stopping criteria tolerance; set to 1e-10 for default
 * @param alpha_ls             for family != 0, line search tuning parameter
 * @param gamma_ls             for family != 0, line search tuning parameter
 * @param max_iter_ls          for family != 0, max number of iterations in line search
 * @param max_iter_newton      for family != 0, max number of iterations in inner ADMM
 * @param DktDk                pointer to the inner product of DktDk
 * @param b                    the link function for a given loss
 * @param b1                   first derivative of the link function for a given loss
 * @param b2                   second derivative of the link function for a given loss
 * @param verbose              0/1 flag for printing progress
 * @return void
 * @see tf_admm
 */
void tf_admm_glm (double * x, double * y, double * w, int n, int k,
		  int max_iter, double lam, int * df,
		  double * beta, double * alpha, double * u,
		  double * obj, int * iter,
		  double rho, double obj_tol, double alpha_ls, double gamma_ls,
		  int max_iter_ls, int max_iter_newton,
		  cs * DktDk, func_RtoR b, func_RtoR b1, func_RtoR b2, int verbose)
{
  double * dir; /* line search direction */
  double * yt;  /* working response: ytilde */
  double * H;   /* weighted Hessian */
  double * obj_admm;

  int i;
  int d;
  int it;
  int iter_admm;
  int * iter_ls;
  double * Db;
  double * Dd;
  double t; /* stepsize */

  dir  = (double*) malloc(n*sizeof(double)); /* line search direction */
  yt = (double*) malloc(n*sizeof(double));   /* working response: ytilde */
  H  = (double*) malloc(n*sizeof(double));   /* weighted Hessian */

  /* Buffers for line search */
  Db      = (double *) malloc(n*sizeof(double));
  Dd      = (double *) malloc(n*sizeof(double));
  iter_ls = (int *)    malloc(sizeof(int));

  obj_admm = (double*) malloc(max_iter*sizeof(double));

  if (verbose) printf("\nlambda=%0.3e\n",lam);
  if (verbose) printf("Iteration\tObjective\tADMM iters\n");

  /* One prox Newton step per iteration */
  for (it=0; it < max_iter_newton; it++)
  {
    /* Define weighted Hessian, and working response */
    for (i=0; i<n; i++)
    {
      H[i] = w[i] * b2(beta[i]);
      if (fabs(H[i])>WEIGHT_SMALL) yt[i] = beta[i] + (y[i]-b1(beta[i]))/H[i];
      else yt[i] = beta[i] + (y[i]-b1(beta[i]));
    }

    /* Prox Newton step */
    iter_admm = 0;
    tf_admm_gauss(x, yt, H, n, k, max_iter, lam, df, dir, alpha, u,
		  obj_admm, &iter_admm, rho, obj_tol, DktDk, 0);

    /* Line search */
    for (i=0; i<n; i++) dir[i] = dir[i] - beta[i];
    t = line_search(x, y, w, n, k, lam, b, b1, beta, dir, alpha_ls, gamma_ls,
		    max_iter_ls, iter_ls, Db, Dd);
    for (i=0; i<n; i++) beta[i] = beta[i] + t * dir[i];

    /* Compute objective */
    obj[it] = tf_obj_glm(x, y, w, n, k, lam, b, beta, yt);
    if (verbose) printf("\t%i\t%0.3e\t%i\n",it+1,obj[it],iter_admm);

    /* Stop if relative difference of objective values < obj_tol */
    if (it > 0 && (fabs(obj[it] - obj[it-1]) < fabs(obj[it-1]) * obj_tol)) break;
  }

  *iter = it;

  /* Compute final df value, based on alpha */
  d = k+1;
  for (i=0; i<n-k; i++) if (alpha[i] != alpha[i+1]) d += 1;
  *df = d;

  /* Free everything */
  free(d);
  free(yt);
  free(H);
  free(iter_ls);
  free(Db);
  free(Dd);
  free(obj_admm);
}

