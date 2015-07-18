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
 * @date 2015-07-14
 * @brief Main functions for fitting fused lasso over a lattice.
 *
 * Here.
 */

#ifndef LATTICE_H
#define LATTICE_H

#include "utils.h"
#include "cs.h"

#define MIN_ITER 10

void do_lattice (double *y, double *w, int n, int m, int p, int max_iter,
                  double lambda, double rho, double eps,
                  int verbose, int naflag,
                  double *beta0,
                  double *beta1, double *beta2, double *beta3,
                  double *thisy1, double *thisy2,  double *thisy3,
                  double *u1, double *u2, double *u3, double *u4,
                  cs *E, double *c, int d,
                  double *buff, double *abuff, int lattice_type);

void do_dp_rows    (double *y, double *buff,
                      double *abuff, double *ans, int n, int m,
                      double lambda);
void do_dp_cols    (double *y, double *buff,
                      double *abuff, double *ans, int n, int m,
                      double lambda);
void do_dp_hexs    (double *y, double *buff,
                      double *abuff, double *ans, int n, int m,
                      double lambda);
void do_dp_cols_na (double *y, double *buff,
                      double *abuff, double *ans, int n, int m,
                      double lambda);
void do_dp_rows_na (double *y, double *buff,
                      double *abuff, double *ans, int n, int m,
                      double lambda);
void do_dp_hexs_na (double *y, double *buff,
                      double *abuff, double *ans, int n, int m,
                      double lambda);

#endif
