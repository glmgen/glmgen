#ifndef TF_H
#define TF_H

void tf_dp(int n, double *y, double lam, double *beta);

void tf_admm (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj, int * numiter,
              double rho, double eabs, double erel);

void tf_prime_dual (double * y, double * x, int n, int k, int family, int max_iter,
              int lam_flag, int obj_flag,  double * lambda, int nlambda, double lambda_min_ratio,
              double * beta, double * obj, int * numiter);

#endif