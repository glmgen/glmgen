#include <Rinternals.h>
#include <string.h>
#include <stdio.h>

#include "tf/tf.h"
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
            SEXP sLamFlag, SEXP sObjFlag, SEXP sLambda, SEXP sNlambda, SEXP sLambdaMinRatio,
            SEXP sControl )
{

  // Initialize all of the variables
  int i;
  double * y;
  double * x;
  double * w;
  int n;
  int k;
  int family;
  int method;
  int maxiter;
  int lam_flag;
  int obj_flag;
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
  SEXP sOutput;
  SEXP sOutputNames;

  double rho;
  int vary_rho;
  double obj_tol;
  // Convert input SEXP variables into C style variables
  y = REAL(sY);
  x = REAL(sX);
  w = REAL(sW);
  n = asInteger(sN);
  k = asInteger(sK);
  family = asInteger(sFamily);
  method = asInteger(sMethod);
  maxiter = asInteger(sMaxIter);
  lam_flag = asInteger(sLamFlag);
  obj_flag = asInteger(sObjFlag);
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

  // Switch on the method, and access low-level C functions
  switch(method)
  {
    case TF_ADMM:
      rho = get_control_value(sControl, "rho", 1);
      vary_rho = get_control_value(sControl, "vary_rho", 0);
      obj_tol = get_control_value(sControl, "obj_tol", 1e-12);

      tf_admm(y, x, w, n, k, family, maxiter, lam_flag, obj_flag, lambda,
              nlambda, lambda_min_ratio, beta, obj, iter, status,
              rho, vary_rho, obj_tol);
      break;

    case TF_PRIMALDUAL_IP:
      tf_primal_dual(y, x, w, n, k, family, maxiter, lam_flag, obj_flag, lambda,
              nlambda, lambda_min_ratio, beta, obj);
      break;

    default:
      error("Method code not found.");
  }

  // Create a list for the output
  PROTECT(sOutput = allocVector(VECSXP, 2 + obj_flag));
  PROTECT(sOutputNames = allocVector(STRSXP, 2 + obj_flag));

  // Assing beta, lambda, and (possibly) obj to the list
  SET_VECTOR_ELT(sOutput, 0, sBeta);
  SET_VECTOR_ELT(sOutput, 1, sLambda);
  if(obj_flag) SET_VECTOR_ELT(sOutput, 2, sObj);

  // Attach names as an attribute to the returned SEXP
  SET_STRING_ELT(sOutputNames, 0, mkChar("beta"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("lambda"));
  if(obj_flag) SET_STRING_ELT(sOutputNames, 2, mkChar("obj"));
  setAttrib(sOutput, R_NamesSymbol, sOutputNames);

  // Free the allocated objects for the gc and return the output as a list
  UNPROTECT(7);
  return sOutput;
}
