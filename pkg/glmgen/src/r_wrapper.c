#include <Rinternals.h>
#include <string.h>
#include <stdio.h>

#include "tf.h"
#include "int_codes.h"

double get_control_value(SEXP sControlList, const char * param_name, double param_default)
{
  SEXP sControlListNames;
  double param_output;
  sControlListNames = getAttrib(sControlList, R_NamesSymbol);
  param_output = param_default;
  for (int i = 0; i < length(sControlList); i++)
  {
    if(strcmp(CHAR(STRING_ELT(sControlListNames, i)), param_name) == 0) {
     param_output = REAL(VECTOR_ELT(sControlList, i))[0];
     //break;
    }
  }
  return param_output;
}

SEXP tf_R ( SEXP sY, SEXP sX, SEXP sN, SEXP sK, SEXP sFamily, SEXP sMethod, SEXP sMaxIter,
            SEXP sLamFlag, SEXP sObjFlag, SEXP sLambda, SEXP sNlambda, SEXP sLambdaMinRatio,
            SEXP sControl )
{

  // Initialize all of the variables
  double * y;
  double * x;
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
  int * numiter;

  SEXP sLambdaNew;
  SEXP sBeta;
  SEXP sObj;
  SEXP sNumIter;
  SEXP sOutput;
  SEXP sOutputNames;

  double rho;
  double eabs;
  double erel;

  // Convert input SEXP variables into C style variables
  y = REAL(sY);
  x = REAL(sX);
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
  PROTECT(sObj = allocVector(REALSXP, maxiter));
  obj = REAL(sObj);
  PROTECT(sNumIter = allocVector(INTSXP, 1));
  numiter = INTEGER(sNumIter);
  numiter[0] = 1;

  // Switch on the method, and access low-level C functions
  switch(method) {
    case TF_ADMM:
      rho = get_control_value(sControl, "rho", 7);
      eabs = get_control_value(sControl, "eabs", 1e-8);
      erel = get_control_value(sControl, "erel", 1e-6);
      tf_admm(y, x, n, k, family, maxiter, lam_flag, obj_flag, lambda,
              nlambda, lambda_min_ratio, beta, obj, numiter,
              rho, eabs, erel);
      break;

    case TF_PRIMALDUAL_IP:
      tf_prime_dual(y, x, n, k, family, maxiter, lam_flag, obj_flag, lambda,
              nlambda, lambda_min_ratio, beta, obj, numiter);
      break;

    default:
      error("Method code not found.");
  }

  // Create a list for the output
  PROTECT(sOutput = allocVector(VECSXP, 2 + 2*obj_flag));
  PROTECT(sOutputNames = allocVector(STRSXP, 2 + 2*obj_flag));

  // Assing beta, lambda, and (possibly) obj and numiter to the list
  SET_VECTOR_ELT(sOutput, 0, sBeta);
  SET_VECTOR_ELT(sOutput, 1, sLambda);
  if(obj_flag) SET_VECTOR_ELT(sOutput, 2, sObj);
  if(obj_flag) SET_VECTOR_ELT(sOutput, 3, sNumIter);

  // Attach names as an attribute to the returned SEXP
  SET_STRING_ELT(sOutputNames, 0, mkChar("beta"));
  SET_STRING_ELT(sOutputNames, 1, mkChar("lambda"));
  if(obj_flag) SET_STRING_ELT(sOutputNames, 2, mkChar("obj"));
  if(obj_flag) SET_STRING_ELT(sOutputNames, 3, mkChar("numiter"));
  setAttrib(sOutput, R_NamesSymbol, sOutputNames);

  // Free the allocated objects for the gc and return the output as a list
  UNPROTECT(6);
  return sOutput;
}
