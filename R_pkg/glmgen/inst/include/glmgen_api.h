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
 * @file glmgen_api.c
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Primary API entry points for user code.
 *
 * Here.
 */

#ifndef GLMGEN_API_H
#define GLMGEN_API_H

 #include "cs.h"
 #include "tf.h"
 #include "utils.h"

/* user-level c functions here */
double * tf_admm_default (double * y, int n);

void tf_admm (double * y, double * x, double * w, int n, int k, int family,
              int max_iter, int lam_flag, double * lambda,
              int nlambda, double lambda_min_ratio, double * beta,
              double * obj, int * iter, int * status, double rho,
              double obj_tol, double alpha_ls, double gamma_ls,
              int max_iter_ls, int max_inner_iter, int verbose);

void tf_predict(double * beta, double * x, int n, int k, int family,
    double * x0, int n0, double * pred, double zero_tol);

void thin( double* x, double* y, double* w, int n, int k,
  double** xt, double** yt, double** wt, int* nt_ptr, double x_cond);

#endif