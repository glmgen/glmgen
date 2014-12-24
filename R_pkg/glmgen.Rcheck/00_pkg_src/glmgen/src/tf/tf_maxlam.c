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
 * @file tf_maxlam.c
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Main calling function for fitting trendfiltering model.
 *
 * Here.
 */

 #include "tf.h"

double tf_maxlam (int len, double * y, gqr * Dt_qr, double * w)
{
  /* This is exact for the Gaussian case, but only approximate
     for logistic or Poisson losses, which is OK, since we are
     only tasked with finding an interesting range for lambdas. */

  int i;
  int tmp; /*n-k-1*/
  double maxlam;
  double * y_work;

  y_work = (double *) malloc(len * sizeof(double));

  for(i = 0; i < len; i++) y_work[i] = sqrt(w[i]) * y[i];

  glmgen_qrsol(Dt_qr, y_work);

  tmp = Dt_qr->n;
  maxlam = 0;
  for(i = 0; i < tmp; i++) maxlam = MAX(maxlam, fabs(y_work[i]));

  free(y_work);
  return maxlam;
}
