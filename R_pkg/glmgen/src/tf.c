#include <math.h>
#include <stdlib.h>

#include "tf.h"
#include "matcomp.h"
#include "utils.h"
#include "int_codes.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

// Need to contruct these (Taylor)
void tf_admm (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda,
              double lambda_min_ratio, double * beta, double * obj,
              double rho, double obj_tol);

void tf_primal_dual (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj);

void tf_admm_gauss (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR);

void tf_admm_logistic (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR);

void tf_admm_pois (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       csn * sparseQR);

void tf_dp (int n, double *y, double lam, double *beta);

cs * tf_calc_Dk (int n, int k, const double * x);

double ts_maxlam (int len, double * Dy, gqr * DDt_qr, int family);

void tf_Dx(double *x, int n, int k,double *a, double *b);

void tf_Dtx(double *x, int n, int k, double *a, double *b);