#ifndef TF_H
#define TF_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "cs.h"
#include "utils.h"
#include "int_codes.h"

#define WEIGHT_SMALL DBL_EPSILON
#define ADMM_MAX_ITER 250

/* Main glmgen api functions */
void tf_admm (double * y, double * x, double * w, int n, int k, int family, int max_iter,
              int lam_flag, double * lambda, int nlambda,
              double lambda_min_ratio, double * beta, double * obj, int * iter,
              int * status, double rho, double obj_tol);

void tf_primal_dual (double * y, double * x, double * w, int n, int k, int family, int max_iter,
              int lam_flag, double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj);

/* Helper functions for cases of the admm and primal dual algorithms */
void tf_dp (int n, double *y, double lam, double *beta);
void tf_dp_weight (int n, double *y, double *w, double lam, double *beta);
void tf_admm_gauss (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, cs * DktDk);
void tf_admm_logistic (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, cs * DktDk);
void tf_admm_poisson (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, cs * DktDk);
       
typedef double (*func_RtoR)(double);
void tf_admm_glm (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, cs * DktDk,
       func_RtoR b, func_RtoR b1, func_RtoR b2);

/* Functions to predict */       
void tf_predict(double * beta, double * x, int n, int k, int family,
	double * x0, int n0, double * pred,
  double zero_tol); 
void tf_predict_gauss(double * beta, double * x, int n, int k,
	    double * x0, int n0, double * pred,
	    double zero_tol);


/* Low-level utility functions for trendfiltering */
cs * tf_calc_dk (int n, int k, const double * x);
cs * tf_calc_dktil (int n, int k, const double * x);
double tf_maxlam (int len, double * y, gqr * Dt_qr, double * w, int family);
void tf_dx(double *x, int n, int k,double *a, double *b);
void tf_dtx(double *x, int n, int k, double *a, double *b);
void tf_dxtil(double *x, int n, int k,double *a, double *b);
void tf_dtxtil(double *x, int n, int k, double *a, double *b);


double tf_line_search(double * y, double * x, double * w, int n, int k, double lam, 
    func_RtoR b, func_RtoR b1, 
    double * beta, double * d, 
    double alpha, double gamma, int max_iter,
    int * iter, double * Db, double * Dd);
#endif
