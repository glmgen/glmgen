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
 * @brief Main functions for fitting trend filtering models.
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

/* ADMM functions */
double * tf_admm_default(double * y, int n);

void tf_admm(double * x, double * y, double * w, int n, int k, int family,
	     int max_iter, int lam_flag, double * lambda, 
	     int nlambda, double lambda_min_ratio, int tridiag, int * df,
	     double * beta, double * obj, int * iter, int * status, 
	     double rho, double obj_tol, double obj_tol_newton, double alpha_ls,
	     double gamma_ls, int max_iter_ls, int max_iter_newton, int verbose);

void tf_admm_gauss(double * x, double * y, double * w, int n, int k,
		   int max_iter, double lam, int * df,
		   double * beta, double * alpha, double * u,
		   double * obj, int * iter,
		   double rho, double obj_tol, cs * DktDk, int verbose);

void tf_admm_gauss_tri (double * x, double * y, double * w, int n, int k,
    int max_iter, double lam, int * df,
    double * beta, double * alpha, double * u,
    double * obj, int * iter,
    double rho, double obj_tol, int verbose);

void tf_admm_glm(double * x, double * y, double * w, int n, int k,
		 int max_iter, double lam, int tridiag, int * df,
		 double * beta, double * alpha, double * u,
		 double * obj, int * iter,
		 double rho, double obj_tol, double obj_tol_newton, double alpha_ls,
		 double gamma_ls, int max_iter_ls, int max_iter_newton,
		 cs * DktDk, func_RtoR b, func_RtoR b1, func_RtoR b2, int verbose);

/* Dynamic programming routines */
void tf_dp (int n, double *y, double lam, double *beta);
void tf_dp_weight (int n, double *y, double *w, double lam, double *beta);

/* Prediction functions */
void tf_predict(double * x, double * beta,  int n, int k, int family,
		double * x0, int n0, double * pred, double zero_tol);
void tf_predict_gauss(double * x, double * beta, int n, int k,
		      double * x0, int n0, double * pred, double zero_tol);
void poly_coefs(double *x, double *beta, int k, double *phi);

/* Calculate maximum lambda for a trend filtering problem */
double tf_maxlam(int n, double * y, gqr * Dt_qr, double * w);

/* Low-level utility functions for trend filtering D matrix */
cs * tf_calc_dk (int n, int k, const double * x);
cs * tf_calc_dktil (int n, int k, const double * x);
void tf_dx (double *x, int n, int k,double *a, double *b);
void tf_dtx (double *x, int n, int k, double *a, double *b);
void tf_dxtil (double *x, int n, int k,double *a, double *b);
void tf_dtxtil (double *x, int n, int k, double *a, double *b);
void tf_dx1(double *x, int n, int j, double *a, double *b);
void tf_dtx1(double *x, int n, int j, double *a, double *b);

/* Compute trend filtering objective */
double tf_obj(double *x, double *y, double *w, int n, int k, double lambda,
	      int family, double *beta, double *buf);
double tf_obj_gauss(double *x, double *y, double *w, int n, int k, double lambda,
		    double *beta, double *buf);
double tf_obj_glm(double *x, double *y, double *w, int n, int k, double lambda,
		  func_RtoR b, double *beta, double *buf);

#endif
