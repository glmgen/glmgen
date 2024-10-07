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
#include <Rdefines.h>
#include <string.h>
#include <stdio.h>

#include "tf.h"
#include "utils.h"
#include "lattice.h"

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
    Rf_error("%s", outputString);
  }

  return param_output;
}

cs * dgTMatrix_to_cs (SEXP input)
{
  int *i;
  int *j;
  int k;
  int m;
  int n;
  int d;
  double *x;
  cs *T;
  cs *U;

  i = INTEGER(GET_SLOT(input,Rf_install("i")));
  j = INTEGER(GET_SLOT(input,Rf_install("j")));
  m = INTEGER(GET_SLOT(input,Rf_install("Dim")))[0];
  n = INTEGER(GET_SLOT(input,Rf_install("Dim")))[1];
  d = LENGTH(GET_SLOT(input,Rf_install("i")));
  x = REAL(GET_SLOT(input,Rf_install("x")));

  T = cs_spalloc (m, n, d, 1, 1);
  for (k = 0; k < d; k++)
  {
    cs_entry (T, i[k], j[k], x[k]);
  }
  U = cs_compress(T);
  cs_spfree(T);
  return U;
}

SEXP thin_R (SEXP sX, SEXP sY, SEXP sW, SEXP sN, SEXP sK, SEXP sControl)
{
  /* Initialize all of the variables */
  int i;
  double * x;
  double * y;
  double * w;
  int n;
  int k;
  double * xt;
  double * yt;
  double * wt;
  int nt;
  double x_tol;

  /* Initalize the output */
  SEXP sOutput;
  SEXP sOutputNames;
  SEXP sOutputX;
  SEXP sOutputY;
  SEXP sOutputW;
  SEXP sOutputN;
  double * outputX;
  double * outputY;
  double * outputW;
  int * outputN;

  /* Grab condition number from the control list */
  x_tol = get_control_value(sControl, "x_tol");

  /* Convert input SEXP variables into C style variables */
  x = REAL(sX);
  y = REAL(sY);
  w = REAL(sW);
  n = asInteger(sN);
  k = asInteger(sK);

  xt = yt = wt = NULL;
  thin(x,y,w,n,k,&xt,&yt,&wt,&nt,x_tol);

  if( xt != NULL )
  {
    x = xt;
    y = yt;
    w = wt;
    n = nt;
  }

  /* Construct the output */
  PROTECT(sOutputX = allocVector(REALSXP, n));
  PROTECT(sOutputY = allocVector(REALSXP, n));
  PROTECT(sOutputW = allocVector(REALSXP, n));
  PROTECT(sOutputN = allocVector(INTSXP, 1));
  outputX = REAL(sOutputX);
  outputY = REAL(sOutputY);
  outputW = REAL(sOutputW);
  outputN = INTEGER(sOutputN);
  memcpy(outputX, x, sizeof(double) * n);
  memcpy(outputY, y, sizeof(double) * n);
  memcpy(outputW, w, sizeof(double) * n);
  memcpy(outputN, &n, sizeof(int));

  PROTECT(sOutput = allocVector(VECSXP, 4));
  PROTECT(sOutputNames = allocVector(STRSXP, 4));
  SET_VECTOR_ELT(sOutput, 0, sOutputX);
  SET_VECTOR_ELT(sOutput, 1, sOutputY);
  SET_VECTOR_ELT(sOutput, 2, sOutputW);
  SET_VECTOR_ELT(sOutput, 3, sOutputN);
  SET_STRING_ELT(sOutputNames, 0, mkChar("x"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("y"));
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

SEXP tf_R ( SEXP sX, SEXP sY, SEXP sW, SEXP sN, SEXP sK, SEXP sFamily, SEXP sMethod,
    SEXP sBeta0, SEXP sLamFlag, SEXP sLambda, SEXP sNlambda, SEXP sLambdaMinRatio,
    SEXP sVerbose, SEXP sControl )
{

  /* Initialize all of the variables */
  int i;
  double * x;
  double * y;
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
  int * df;
  double * beta;
  double * obj;
  int * iter;
  int * status;
  int verbose;
  double * beta0;

  SEXP sLambdaNew;
  SEXP sDf;
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
  double obj_tol_newton;
  double alpha_ls;
  double gamma_ls;
  int max_iter_ls;
  int max_iter_newton;
	int max_iter_outer;
	int tridiag;			

  /* Convert input SEXP variables into C style variables */
  x = REAL(sX);
  y = REAL(sY);
  w = REAL(sW);
  n = asInteger(sN);
  k = asInteger(sK);
  verbose = asInteger(sVerbose);

  family = asInteger(sFamily);
  method = asInteger(sMethod);
  max_iter = get_control_value(sControl, "max_iter");
  max_iter_newton = get_control_value(sControl, "max_iter_newton");
	max_iter_outer = (family == FAMILY_GAUSSIAN) ? max_iter : max_iter_newton;
  lam_flag = asInteger(sLamFlag);
  PROTECT(sLambdaNew = duplicate(sLambda));
  lambda = REAL(sLambda);
  nlambda = asInteger(sNlambda);
  lambda_min_ratio = asReal(sLambdaMinRatio);
  beta0 = isNull(sBeta0) ? NULL : REAL(sBeta0);
  
  PROTECT(sDf = allocVector(INTSXP, nlambda));
  df = INTEGER(sDf);
  for(i = 0; i < nlambda; i++) df[i] = 0;
  PROTECT(sBeta = allocMatrix(REALSXP, n, nlambda));
  beta = REAL(sBeta);

  PROTECT(sObj = allocMatrix(REALSXP, max_iter_outer+1, nlambda));
  obj = REAL(sObj);
  for(i = 0; i < (max_iter_outer+1) * nlambda; i++) obj[i] = 0;
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
      obj_tol_newton = get_control_value(sControl, "obj_tol_newton");
      alpha_ls = get_control_value(sControl, "alpha_ls");
      gamma_ls = get_control_value(sControl, "gamma_ls");
      max_iter_ls = get_control_value(sControl, "max_iter_ls");
      tridiag = get_control_value(sControl, "tridiag");
			
      tf_admm(x, y, w, n, k, family, max_iter, beta0, lam_flag, lambda,
          nlambda, lambda_min_ratio, tridiag, df, beta, obj, iter, status,
          rho, obj_tol, obj_tol_newton, alpha_ls, gamma_ls, max_iter_ls,
          max_iter_newton, verbose);
      break;

    default:
      error("Method code not found.");
      break;
  }

  /* Create a list for the output */
  PROTECT(sOutput = allocVector(VECSXP, 9));
  PROTECT(sOutputNames = allocVector(STRSXP, 9));

  /* Adding beta, lambda, and obj to the list */
  SET_VECTOR_ELT(sOutput, 0, sBeta);
  SET_VECTOR_ELT(sOutput, 1, sLambda);
  SET_VECTOR_ELT(sOutput, 2, sDf);
  SET_VECTOR_ELT(sOutput, 3, sObj);
  SET_VECTOR_ELT(sOutput, 4, sIter);
  SET_VECTOR_ELT(sOutput, 5, sStatus);
  SET_VECTOR_ELT(sOutput, 6, sXt);
  SET_VECTOR_ELT(sOutput, 7, sYt);
  SET_VECTOR_ELT(sOutput, 8, sWt);


  /* Attach names as an attribute to the returned SEXP */
  SET_STRING_ELT(sOutputNames, 0, mkChar("beta"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("lambda"));
  SET_STRING_ELT(sOutputNames, 2, mkChar("df"));
  SET_STRING_ELT(sOutputNames, 3, mkChar("obj"));
  SET_STRING_ELT(sOutputNames, 4, mkChar("iter"));
  SET_STRING_ELT(sOutputNames, 5, mkChar("status"));
  SET_STRING_ELT(sOutputNames, 6, mkChar("x"));
  SET_STRING_ELT(sOutputNames, 7, mkChar("y"));
  SET_STRING_ELT(sOutputNames, 8, mkChar("w"));
  setAttrib(sOutput, R_NamesSymbol, sOutputNames);

  /* Free the allocated objects for the gc and return the output as a list */
  UNPROTECT(11);
  return sOutput;
}

SEXP tf_predict_R (SEXP sX, SEXP sBeta, SEXP sN, SEXP sK, SEXP sX0, SEXP sN0,
    SEXP sNLambda, SEXP sFamily, SEXP sZeroTol)
{
  /* Initialize all of the variables */
  int i;
  double * x;
  double * beta;
  double * x0;
  int n;
  int n0;
  int k;
  int family;
  int nlambda;
  double zero_tol;
  double * pred;

  x = REAL(sX);
  beta = REAL(sBeta);
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
    tf_predict(x, beta+i*n, n, k, family, x0, n0, pred + n0*i, zero_tol);
  }

  /* Free the allocated objects for the gc and return the output */
  UNPROTECT(1);
  return sPred;
}

SEXP lattice_R (SEXP sY, SEXP sW, SEXP sWedge, SEXP sLambda, SEXP sRho,
    SEXP sEps, SEXP sMaxiter, SEXP sVerbose, SEXP sNaflag,
    SEXP sLatticeType, SEXP sMethodType,
    SEXP sE, SEXP sC, SEXP sBeta0)
{
  int k;
  int n;
  int m;
  int p;
  int d;
  int N;
  int edge_length;
  int max_iter;
  int verbose;
  int naflag;
  int lattice_type;
  int method_type;
  double thisLam;
  double rho;
  double eps;
  double *y;
  double *w;
  double *ew;
  double *lambda;
  double *output;
  double *buff;
  double *abuff;
  double *wbuff;
  double *thisy1;
  double *thisy2;
  double *thisy3;
  double *thisy4;
  double *beta0;
  double *beta1;
  double *beta2;
  double *beta3;
  double *u1;
  double *u2;
  double *u3;
  double *u4;
  double *c;
  cs *E;
  SEXP sOutput;

  /* Load inputs into native c types*/
  n = INTEGER(GET_DIM(sY))[0];
  m = INTEGER(GET_DIM(sY))[1];
  max_iter = INTEGER(sMaxiter)[0];
  verbose = INTEGER(sVerbose)[0];
  naflag = INTEGER(sNaflag)[0];
  lattice_type = INTEGER(sLatticeType)[0];
  method_type = INTEGER(sMethodType)[0];
  if (lattice_type == LATTICE_3D_GRID)
  {
    p = INTEGER(GET_DIM(sY))[2];
  } else {
    p = 1;
  };
  thisLam = REAL(sLambda)[0];
  rho = REAL(sRho)[0];
  eps = REAL(sEps)[0];
  y = REAL(sY);
  w = REAL(sW);
  ew = (TYPEOF(sWedge) == NILSXP) ? NULL : REAL(sWedge);
  N = n * m * p;

  /* Parse the (optional) constraint Eb=c */
  E = dgTMatrix_to_cs(sE);
  c = REAL(sC);
  d = E->m;

  /* Determine number of edges */
  edge_length = (n-1)*m*p + n*(m-1)*p;
  switch(lattice_type)
  {
    case LATTICE_2D_GRID:
      break;

    case LATTICE_HEX_GRID:
      edge_length += (n-1)*(m-1)*p;
      break;

    case LATTICE_3D_GRID:
      edge_length += n*m*(p-1);
      break;
  }

  /* Allocate output and working buffers */
  sOutput = PROTECT(allocVector(REALSXP, N));
  output  = REAL(sOutput);
  buff    = REAL(PROTECT(allocVector(REALSXP, N)));
  abuff   = REAL(PROTECT(allocVector(REALSXP, MAX(MAX(n,m), p))));
  wbuff   = REAL(PROTECT(allocVector(REALSXP, MAX(MAX(n,m), p))));
  thisy1  = REAL(PROTECT(allocVector(REALSXP, N)));
  thisy2  = REAL(PROTECT(allocVector(REALSXP, N)));
  thisy3  = REAL(PROTECT(allocVector(REALSXP, N)));
  thisy4  = REAL(PROTECT(allocVector(REALSXP, N)));
  beta0   = REAL(PROTECT(allocVector(REALSXP, N)));
  beta1   = REAL(PROTECT(allocVector(REALSXP, N)));
  beta2   = REAL(PROTECT(allocVector(REALSXP, N)));
  beta3   = REAL(PROTECT(allocVector(REALSXP, N)));
  u1      = REAL(PROTECT(allocVector(REALSXP, N)));
  u2      = REAL(PROTECT(allocVector(REALSXP, N)));
  u3      = REAL(PROTECT(allocVector(REALSXP, N)));
  u4      = REAL(PROTECT(allocVector(REALSXP, d)));
  lambda  = REAL(PROTECT(allocVector(REALSXP, edge_length+1)));

  for (k = 0; k < N; k++)
  {
    beta0[k] = beta1[k] = beta2[k] = beta3[k] = REAL(sBeta0)[k];
    u1[k] = u2[k] = u3[k] = 0;
  }
  for (k = 0; k < d; k++)
  {
    u4[k] = 0;
  }
  for (k = 0; k < edge_length; k++)
  {
    lambda[k] = (ew == NULL) ? thisLam : thisLam*ew[k];
  }
  lambda[edge_length] = 1;

  do_lattice(y, w, n, m, p, max_iter, lambda, rho, eps,
      verbose, naflag,
      beta0, beta1, beta2, beta3,
      thisy1, thisy2, thisy3, thisy4,
      u1, u2, u3, u4,
      E, c, d,
      buff, abuff, wbuff,
      lattice_type, method_type);

  memcpy(output, beta0, sizeof(double) * n * m * p);

  UNPROTECT(17);
  return sOutput;
}


SEXP graph_fused_R (SEXP sY, SEXP sW, SEXP sEdge, SEXP sWedge,
    SEXP sEdgeLen,
    SEXP sLambda, SEXP sRho, SEXP sEps, SEXP sMaxiter,
    SEXP sVerbose, SEXP sMethodType,
    SEXP sE, SEXP sC, SEXP sBeta0)
{
  int i;
  int k;
  int n;
  int num_chains;
  int num_edge_index;
  int d;
  int max_iter;
  int verbose;
  int method_type;
  int *elen;
  int *e;
  int  *ebuff;
  double thisLam;
  double rho;
  double eps;
  double *y;
  double *w;
  double *ew;
  double *lambda;
  double *output;
  double *buff;
  double *abuff;
  double *wbuff;
  double *thisY;
  double *beta0;
  double *B;
  double *U;
  double *u4;
  double *c;
  cs *E;
  SEXP sOutput;

  /* Load inputs into native c types*/
  n = LENGTH(sY);
  num_chains = LENGTH(sEdgeLen);
  num_edge_index = LENGTH(sEdge);
  max_iter = INTEGER(sMaxiter)[0];
  verbose = INTEGER(sVerbose)[0];
  method_type = INTEGER(sMethodType)[0];
  thisLam = REAL(sLambda)[0];
  rho = REAL(sRho)[0];
  eps = REAL(sEps)[0];
  elen = INTEGER(sEdgeLen);
  y = REAL(sY);
  w = REAL(sW);
  e = INTEGER(sEdge);
  ew = REAL(sWedge);

  /* Parse the (optional) constraint Eb=c */
  E = dgTMatrix_to_cs(sE);
  c = REAL(sC);
  d = E->m;

  /* Allocate output and working buffers */
  sOutput = PROTECT(allocVector(REALSXP, n));
  output  = REAL(sOutput);
  buff    = REAL(PROTECT(allocVector(REALSXP, n)));
  abuff   = REAL(PROTECT(allocVector(REALSXP, n)));
  wbuff   = REAL(PROTECT(allocVector(REALSXP, n)));
  ebuff   = INTEGER(PROTECT(allocVector(INTSXP, n)));
  thisY   = REAL(PROTECT(allocVector(REALSXP, n*(num_chains+1))));
  beta0   = REAL(PROTECT(allocVector(REALSXP, n)));
  B       = REAL(PROTECT(allocVector(REALSXP, n*num_chains)));
  U       = REAL(PROTECT(allocVector(REALSXP, n*num_chains)));
  u4      = REAL(PROTECT(allocVector(REALSXP, d)));
  lambda  = REAL(PROTECT(allocVector(REALSXP, num_edge_index+1)));

  for (k = 0; k < n; k++)
  {
    beta0[k] = REAL(sBeta0)[k];
    for (i = 0; i < num_chains; i++)
    {
      B[k + i*n] = REAL(sBeta0)[k];
      U[k + i*n] = 0;
    }
  }
  for (k = 0; k < d; k++)
  {
    u4[k] = 0;
  }
  for (k = 0; k < num_edge_index; k++)
  {
    lambda[k] = thisLam*ew[k];
  }
  lambda[num_edge_index] = 1;

  do_fused_graph(y, w, e, elen, n, num_chains, max_iter,
      lambda, rho, eps,
      verbose, beta0, B, thisY, U, u4, E, c, d,
      buff, abuff, wbuff, ebuff,
      method_type);

  memcpy(output, beta0, sizeof(double) * n);

  UNPROTECT(11);
  return sOutput;
}

SEXP matMultiply_R (SEXP sB, SEXP sK, SEXP sX, SEXP sMatrixCode)
{
  double *b;
  double *x;
  int k;
  int matcode;
  int n;

  SEXP sOutput;
  double *output;

  b = REAL(sB);
  k = INTEGER(sK)[0];
  x = REAL(sX);
  matcode = INTEGER(sMatrixCode)[0];
  n = LENGTH(sB);
  sOutput = PROTECT(allocVector(REALSXP, LENGTH(sB)));
  output = REAL(sOutput);

  switch(matcode)
  {
    case 0:
      tf_dx(x, n, k, b, output);
      break;

    default:
      error("Method code not found.");
      break;
  }

  UNPROTECT(1);
  return sOutput;
}






