/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Veeranjaneyulu Sadhanala,           *
 *                       Ryan Tibshirani                                    *
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

#include <Rinternals.h>
#include <string.h>
#include <stdio.h>

#include "tf.h"
#include "utils.h"

double get_control_value(SEXP sControlList, const char * param_name)
{
  int i;
  int okay_flag;
  double param_output;
  char * outputString;
  SEXP sControlListNames;

  sControlListNames = getAttrib(sControlList, R_NamesSymbol);
  okay_flag = 0;
  param_output = 0;
  for (i = 0; i < length(sControlList); i++)
  {
    if (strcmp(CHAR(STRING_ELT(sControlListNames, i)), param_name) == 0)
    {
     param_output = REAL(VECTOR_ELT(sControlList, i))[0];
     okay_flag = 1;
     break;
    }
  }

  if (okay_flag == 0)
  {
    outputString = (char *) malloc(sizeof(char) * 100);
    snprintf(outputString, sizeof(char) * 100, "Missing required tuning parameter %s", param_name);
    Rf_error(outputString);
  }

  return param_output;
}

SEXP thin_R (SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sControl)
{
  /* Initialize all of the variables */
  int i;
  double * y;
  double * x;
  double * w;
  int n;
  int k;
  double * xt;
  double * yt;
  double * wt;
  int nt;
  double x_cond;

  /* Initalize the output */
  SEXP sOutput;
  SEXP sOutputNames;
  SEXP sOutputY;
  SEXP sOutputX;
  SEXP sOutputW;
  SEXP sOutputN;
  double * outputY;
  double * outputX;
  double * outputW;
  int * outputN;

  /* Grab condition number from the control list */
  x_cond = get_control_value(sControl, "x_cond");

  /* Convert input SEXP variables into C style variables */
  y = REAL(sY);
  x = REAL(sX);
  w = REAL(sW);
  n = asInteger(sN);
  k = asInteger(sK);

  xt = yt = wt = NULL;
  thin(x,y,w,n,k,&xt,&yt,&wt,&nt,x_cond);

  if( xt != NULL )
  {
    x = xt;
    y = yt;
    w = wt;
    n = nt;
  }

  /* Construct the output */
  PROTECT(sOutputY = allocVector(REALSXP, n));
  PROTECT(sOutputX = allocVector(REALSXP, n));
  PROTECT(sOutputW = allocVector(REALSXP, n));
  PROTECT(sOutputN = allocVector(INTSXP, 1));
  outputY = REAL(sOutputY);
  outputX = REAL(sOutputX);
  outputW = REAL(sOutputW);
  outputN = INTEGER(sOutputN);
  memcpy(outputY, y, sizeof(double) * n);
  memcpy(outputX, x, sizeof(double) * n);
  memcpy(outputW, w, sizeof(double) * n);
  memcpy(outputN, &n, sizeof(int));

  PROTECT(sOutput = allocVector(VECSXP, 4));
  PROTECT(sOutputNames = allocVector(STRSXP, 4));
  SET_VECTOR_ELT(sOutput, 0, sOutputY);
  SET_VECTOR_ELT(sOutput, 1, sOutputX);
  SET_VECTOR_ELT(sOutput, 2, sOutputW);
  SET_VECTOR_ELT(sOutput, 3, sOutputN);
  SET_STRING_ELT(sOutputNames, 0, mkChar("y"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("x"));
  SET_STRING_ELT(sOutputNames, 2, mkChar("w"));
  SET_STRING_ELT(sOutputNames, 3, mkChar("n"));
  setAttrib(sOutput, R_NamesSymbol, sOutputNames);

  /* Free temporary thinning variables */
  if(xt != NULL) free(xt);
  if(yt != NULL) free(yt);
  if(wt != NULL) free(wt);

  // Free the allocated objects for the gc and return the output as a list
  UNPROTECT(6);
  return sOutput;
}

SEXP tf_R ( SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sFamily, SEXP sMethod,
            SEXP sLamFlag, SEXP sLambda, SEXP sNlambda, SEXP sLambdaMinRatio,
            SEXP sVerbose, SEXP sControl )
{

  /* Initialize all of the variables */
  int i;
  double * y;
  double * x;
  double * w;
  double * sxt;
  double * syt;
  double * swt;

  int n;
  int k;
  int nt;
  int family;
  int method;
  int max_iter;
  int lam_flag;
  double * lambda;
  int nlambda;
  double lambda_min_ratio;
  double * beta;
  double * obj;
  int * iter;
  int * status;
  int verbose;

  SEXP sLambdaNew;
  SEXP sBeta;
  SEXP sObj;
  SEXP sIter;
  SEXP sStatus;
  SEXP sXt;
  SEXP sYt;
  SEXP sWt;
  SEXP sOutput;
  SEXP sOutputNames;

  double rho;
  double obj_tol;
  double alpha_ls;
  double gamma_ls;
  int max_iter_ls;
  int max_iter_newton;

  /* Convert input SEXP variables into C style variables */
  y = REAL(sY);
  x = REAL(sX);
  w = REAL(sW);
  n = asInteger(sN);
  k = asInteger(sK);
  verbose = asInteger(sVerbose);

  family = asInteger(sFamily);
  method = asInteger(sMethod);
  max_iter = get_control_value(sControl, "max_iter");
  lam_flag = asInteger(sLamFlag);
  PROTECT(sLambdaNew = duplicate(sLambda));
  lambda = REAL(sLambda);
  nlambda = asInteger(sNlambda);
  lambda_min_ratio = asReal(sLambdaMinRatio);
  PROTECT(sBeta = allocMatrix(REALSXP, n, nlambda));
  beta = REAL(sBeta);
  PROTECT(sObj = allocMatrix(REALSXP, max_iter, nlambda));
  obj = REAL(sObj);
  for(i = 0; i < max_iter * nlambda; i++) obj[i] = 0;
  PROTECT(sIter = allocVector(INTSXP, nlambda));
  iter = INTEGER(sIter);
  for(i = 0; i < nlambda; i++) iter[i] = 0;
  PROTECT(sStatus = allocVector(INTSXP, nlambda));
  status = INTEGER(sStatus);
  for(i = 0; i < nlambda; i++) status[i] = 0;

  PROTECT(sXt = allocVector(REALSXP, n));
  sxt = REAL(sXt);
  for(i = 0; i < n; i++) sxt[i] = x[i];
  PROTECT(sYt = allocVector(REALSXP, n));
  syt = REAL(sYt);
  for(i = 0; i < n; i++) syt[i] = y[i];
  PROTECT(sWt = allocVector(REALSXP, n));
  swt = REAL(sWt);
  for(i = 0; i < n; i++) swt[i] = w[i];

  /* Switch on the method, and access low-level C functions */
  switch(method)
  {
    case TF_ADMM:
      rho = get_control_value(sControl, "rho");
      obj_tol = get_control_value(sControl, "obj_tol");
      alpha_ls = get_control_value(sControl, "alpha_ls");
      gamma_ls = get_control_value(sControl, "gamma_ls");
      max_iter_ls = get_control_value(sControl, "max_iter_ls");
      max_iter_newton = get_control_value(sControl, "max_iter_newton");

      tf_admm(y, x, w, n, k, family, max_iter, lam_flag, lambda,
              nlambda, lambda_min_ratio, beta, obj, iter, status,
              rho, obj_tol, alpha_ls, gamma_ls, max_iter_ls,
              max_iter_newton, verbose);
      break;

    default:
      error("Method code not found.");
      break;
  }

  /* Create a list for the output */
  PROTECT(sOutput = allocVector(VECSXP, 8));
  PROTECT(sOutputNames = allocVector(STRSXP, 8));

  /* Assing beta, lambda, and obj to the list */
  SET_VECTOR_ELT(sOutput, 0, sBeta);
  SET_VECTOR_ELT(sOutput, 1, sLambda);
  SET_VECTOR_ELT(sOutput, 2, sObj);
  SET_VECTOR_ELT(sOutput, 3, sIter);
  SET_VECTOR_ELT(sOutput, 4, sStatus);
  SET_VECTOR_ELT(sOutput, 5, sXt);
  SET_VECTOR_ELT(sOutput, 6, sYt);
  SET_VECTOR_ELT(sOutput, 7, sWt);


  /* Attach names as an attribute to the returned SEXP */
  SET_STRING_ELT(sOutputNames, 0, mkChar("beta"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("lambda"));
  SET_STRING_ELT(sOutputNames, 2, mkChar("obj"));
  SET_STRING_ELT(sOutputNames, 3, mkChar("iter"));
  SET_STRING_ELT(sOutputNames, 4, mkChar("status"));
  SET_STRING_ELT(sOutputNames, 5, mkChar("x"));
  SET_STRING_ELT(sOutputNames, 6, mkChar("y"));
  SET_STRING_ELT(sOutputNames, 7, mkChar("w"));
  setAttrib(sOutput, R_NamesSymbol, sOutputNames);

  /* Free the allocated objects for the gc and return the output as a list */
  UNPROTECT(10);
  return sOutput;
}

SEXP tf_predict_R (SEXP sBeta, SEXP sX, SEXP sN, SEXP sK, SEXP sX0, SEXP sN0,
    SEXP sNLambda, SEXP sFamily, SEXP sZeroTol)
{
  /* Initialize all of the variables */
  int i;
  double * beta;
  double * x;
  double * x0;
  int n;
  int n0;
  int k;
  int family;
  int nlambda;
  double zero_tol;
  double * pred;

  beta = REAL(sBeta);
  x = REAL(sX);
  x0 = REAL(sX0);
  n = asInteger(sN);
  n0 = asInteger(sN0);
  k = asInteger(sK);
  nlambda = asInteger(sNLambda);
  family = asInteger(sFamily);
  zero_tol = asReal(sZeroTol);

  /* Output */
  SEXP sPred;
  PROTECT(sPred = allocVector(REALSXP, n0 * nlambda));
  pred = REAL(sPred);

  for (i = 0; i < nlambda; i++)
  {
    tf_predict(beta, x, n, k, family, x0, n0, pred + n0*i, zero_tol);
  }

  /* Free the allocated objects for the gc and return the output as a list */
  UNPROTECT(1);
  return sPred;
}





