#include <Rinternals.h>
#include <string.h>
#include <stdio.h>

#include "tf/tf.h"
#include "utils/utils.h"
#include "utils/int_codes.h"

double get_control_value(SEXP sControlList, const char * param_name, double param_default)
{
  int i;
  double param_output;
  SEXP sControlListNames;

  sControlListNames = getAttrib(sControlList, R_NamesSymbol);
  param_output = param_default;
  for (i = 0; i < length(sControlList); i++)
  {
    if (strcmp(CHAR(STRING_ELT(sControlListNames, i)), param_name) == 0)
    {
     param_output = REAL(VECTOR_ELT(sControlList, i))[0];
     break;
    }
  }

  return param_output;
}

SEXP tf_R ( SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sFamily, SEXP sMethod, SEXP sMaxIter,
            SEXP sLamFlag, SEXP sLambda, SEXP sNlambda, SEXP sLambdaMinRatio, SEXP sControl )
{

  // Initialize all of the variables
  int i;
  double * y;
  double * x;
  double * w;
  double * xt, * sxt;
  double * yt, * syt;
  double * wt, * swt;

  int n;
  int k;
  int nt;
  int family;
  int method;
  int maxiter;
  int lam_flag;
  double * lambda;
  int nlambda;
  double lambda_min_ratio;
  double * beta;
  double * obj;
  int * iter;
  int * status;
  

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
  double do_thin;
  double x_cond;

  // Convert input SEXP variables into C style variables
  y = REAL(sY);
  x = REAL(sX);
  w = REAL(sW);
  n = asInteger(sN);
  k = asInteger(sK);

  do_thin = get_control_value(sControl, "thinning", 0);
  x_cond = get_control_value(sControl, "x_cond", 0);
  
  xt = yt = wt = NULL;
  if( do_thin > 0 || x_cond > 0 )
  {    
    thin(x,y,w,n,k,&xt,&yt,&wt,&nt,x_cond);
    int thinned = (xt != NULL);
    if( thinned ) 
    {
      x = xt;
      y = yt;
      w = wt;
      n = nt;
    }
  }

  family = asInteger(sFamily);
  method = asInteger(sMethod);
  maxiter = asInteger(sMaxIter);
  lam_flag = asInteger(sLamFlag);
  PROTECT(sLambdaNew = duplicate(sLambda));
  lambda = REAL(sLambda);
  nlambda = asInteger(sNlambda);
  lambda_min_ratio = asReal(sLambdaMinRatio);
  PROTECT(sBeta = allocMatrix(REALSXP, n, nlambda));
  beta = REAL(sBeta);
  PROTECT(sObj = allocMatrix(REALSXP, maxiter, nlambda));
  obj = REAL(sObj);
  for(i = 0; i < maxiter * nlambda; i++) obj[i] = 0;
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

  // Switch on the method, and access low-level C functions
  switch(method)
  {
    case TF_ADMM:
      rho = get_control_value(sControl, "rho", 1);
      obj_tol = get_control_value(sControl, "obj_tol", 1e-10);

      tf_admm(y, x, w, n, k, family, maxiter, lam_flag, lambda,
          nlambda, lambda_min_ratio, beta, obj, iter, status,
          rho, obj_tol);
      break;

    case TF_PRIMALDUAL_IP:
      tf_primal_dual(y, x, w, n, k, family, maxiter, lam_flag, lambda,
          nlambda, lambda_min_ratio, beta, obj);
      break;

    default:
      error("Method code not found.");
  }

  // Create a list for the output
  PROTECT(sOutput = allocVector(VECSXP, 7));
  PROTECT(sOutputNames = allocVector(STRSXP, 7));

  // Assing beta, lambda, and obj to the list
  SET_VECTOR_ELT(sOutput, 0, sBeta);
  SET_VECTOR_ELT(sOutput, 1, sLambda);
  SET_VECTOR_ELT(sOutput, 2, sObj);
  SET_VECTOR_ELT(sOutput, 3, sIter);
  SET_VECTOR_ELT(sOutput, 4, sXt);
  SET_VECTOR_ELT(sOutput, 5, sYt);
  SET_VECTOR_ELT(sOutput, 6, sWt);


  // Attach names as an attribute to the returned SEXP
  SET_STRING_ELT(sOutputNames, 0, mkChar("beta"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("lambda"));
  SET_STRING_ELT(sOutputNames, 2, mkChar("obj"));
  SET_STRING_ELT(sOutputNames, 3, mkChar("iter"));
  SET_STRING_ELT(sOutputNames, 4, mkChar("x"));
  SET_STRING_ELT(sOutputNames, 5, mkChar("y"));
  SET_STRING_ELT(sOutputNames, 6, mkChar("w"));
  setAttrib(sOutput, R_NamesSymbol, sOutputNames);

  if(xt != NULL) free(xt);
  if(yt != NULL) free(yt);
  if(wt != NULL) free(wt);
  // Free the allocated objects for the gc and return the output as a list
  UNPROTECT(10);
  return sOutput;
}
