#ifndef TF_H
#define TF_H
#ifndef _CS_H
#include "cs.h"
#endif

// Primary function for calling the admm algorithm
void tf_admm (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj,
              double rho, double obj_tol);

// Primary function for calling the prime dual algorithm
void tf_primal_dual (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj);

// These compute one step of the admm algorithm for the given family
void tf_admm_gauss(double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR);
void tf_admm_logistic(double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR);
void tf_admm_pois(double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR);

// Functions to calculate various quantities before calling ts_admm_FAMILY
double ts_maxlam(double * y, double * x, int n, int k, double * beta_max);
void tf_calc_dtd(double * x, int n, int k, double * dtd);
void tf_getrho(double * rho, double lambda);
void tf_calc_sparse_qr(int n, int k, double rho, double * dtd, csn * sparseQR);

// Lower level helper functions
void tf_dp(int n, double *y, double lam, double *beta);
void tf_d(double *x, int n, int k,double *a, double *b);
void tf_dt(double *x, int n, int k, double *a, double *b);

#endif
