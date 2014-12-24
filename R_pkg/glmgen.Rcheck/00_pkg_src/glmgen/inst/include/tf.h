/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala   *
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
 * @file tf.h
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Main calling function for fitting trendfiltering model.
 *
 * Here.
 */

#ifndef TF_H
#define TF_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "cs.h"
#include "utils.h"

#define WEIGHT_SMALL DBL_EPSILON

/* Main glmgen api functions */
void tf_admm (double * y, double * x, double * w, int n, int k, int family,
              int max_iter, int lam_flag, double * lambda,
              int nlambda, double lambda_min_ratio, double * beta,
              double * obj, int * iter, int * status, double rho,
              double obj_tol, double alpha_ls, double gamma_ls,
              int max_iter_ls, int max_inner_iter, int verbose);

/* Helper functions for cases of the admm and primal dual algorithms */
void tf_dp (int n, double *y, double lam, double *beta);
void tf_dp_weight (int n, double *y, double *w, double lam, double *beta);
void tf_admm_gauss (double * y, double * x, double * w, int n, int k,
    int max_iter, double lam,
    double * beta, double * alpha, double * u,
    double * obj, int * iter,
    double rho, double obj_tol, cs * DktDk, int verbose);
void tf_admm_glm (double * y, double * x, double * w, int n, int k,
    int max_iter, double lam,
    double * beta, double * alpha, double * u,
    double * obj, int * iter,
    double rho, double obj_tol, double alpha_ls, double gamma_ls,
    int max_iter_ls, int max_iter_admm,
    cs * DktDk,
    func_RtoR b, func_RtoR b1, func_RtoR b2, int verbose);

/* Functions to predict */
void tf_predict(double * beta, double * x, int n, int k, int family,
    double * x0, int n0, double * pred, double zero_tol);
void tf_predict_gauss(double * beta, double * x, int n, int k,
    double * x0, int n0, double * pred, double zero_tol);

/* Low-level utility functions for trendfiltering */
cs * tf_calc_dk (int n, int k, const double * x);
cs * tf_calc_dktil (int n, int k, const double * x);
double tf_maxlam (int len, double * y, gqr * Dt_qr, double * w);
void tf_dx(double *x, int n, int k,double *a, double *b);
void tf_dtx(double *x, int n, int k, double *a, double *b);
void tf_dxtil(double *x, int n, int k,double *a, double *b);
void tf_dtxtil(double *x, int n, int k, double *a, double *b);
void poly_coefs(double *x, int k, double *beta, double *phi);


#endif
