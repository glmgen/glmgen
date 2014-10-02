#ifndef TF_H
#define TF_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cs.h"
#include "tf.h"
#include "utils.h"
#include "int_codes.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* Main glmgen api functions */
void tf_admm (double * y, double * x, double * w, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda,
              double lambda_min_ratio, double * beta, double * obj,
              double rho, double obj_tol);

void tf_primal_dual (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj);

/* Helper functions for cases of the admm and primal dual algorithms */
void tf_dp (int n, double *y, double lam, double *beta);
void tf_admm_gauss (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR);
void tf_admm_logistic (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR);
void tf_admm_pois (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR);
typedef double (*func_RtoR)(double);
void tf_admm_glm (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR,
       func_RtoR f1, func_RtoR f2);

/* Low-level utility functions for trendfiltering */
cs * tf_calc_dk (int n, int k, const double * x);
double ts_maxlam (int len, double * Dy, gqr * DDt_qr, int family);
void tf_dx(double *x, int n, int k,double *a, double *b);
void tf_dtx(double *x, int n, int k, double *a, double *b);

#endif
