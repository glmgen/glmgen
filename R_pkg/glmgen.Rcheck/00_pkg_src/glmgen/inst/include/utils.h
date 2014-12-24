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
 * @file utils.h
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Main calling function for fitting trendfiltering model.
 *
 * Here.
 */

#ifndef UTILS_H
#define UTILS_H

#include "cs.h"

#define TF_ADMM 0
#define TF_PRIMALDUAL_IP 1
#define TF_PROJECTED_NEWTON 2

#define FAMILY_GAUSSIAN 0
#define FAMILY_LOGISTIC 1
#define FAMILY_POISSON 2

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* --- custom structure to hold the out of the qr decomposition ------- */
typedef struct gcs_qr  /* holds numeric and symbol qr together */
{
  csi m ;         /* number of rows */
  csi n ;         /* number of columns */
  css *S ;        /* symbolic qr output */
  csn *N ;        /* numeric qr output */
  double * W ;    /* pre-allocated working space for the solver */
} gqr ;

typedef double (*func_RtoR)(double);

/* Utility functions for solving a linear system with a
   pre-computed sparse qr object (struct gcs_qr / gqr) */
gqr * glmgen_qr (const cs * A);
csi glmgen_qrsol (gqr * B, double * b);
csi glmgen_gqr_free (gqr * A);

cs * scalar_plus_diag (const cs * A, double b, double *D);
void diag_times_sparse (const cs * A, double * w);
double glmgen_factorial(int n);

double glmgen_factorial(int n);
double l1norm(double * x, int n);
int is_nan(double x);
int has_nan(double * x, int n);
int count_nans(double * x, int n);

void thin( double* x, double* y, double* w, int n, int k,
  double** xt, double** yt, double** wt, int* nt_ptr, double x_cond);

double line_search(double * y, double * x, double * w, int n, int k, double lam,
    func_RtoR b, func_RtoR b1,
    double * beta, double * d,
    double alpha, double gamma, int max_iter,
    int * iter, double * Db, double * Dd);

/* glm loss functions */
double logi_b(double x);
double logi_b1(double x);
double logi_b2(double x);
double pois_b(double x);
double pois_b1(double x);
double pois_b2(double x);

#endif
