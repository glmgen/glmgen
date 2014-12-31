#ifndef RWRAPPER_H
#define RWRAPPER_H

double get_control_value(SEXP sControlList, const char * param_name);

SEXP thin_R (SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sControl);
SEXP tf_R (SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sFamily,
    SEXP sMethod, SEXP sLamFlag, SEXP sLambda, SEXP sNlambda,
    SEXP sLambdaMinRatio, SEXP sControl );
SEXP tf_predict_R (SEXP sBeta, SEXP sX, SEXP sN, SEXP sK, SEXP sX0, SEXP sN0,
    SEXP sNLambda, SEXP sFamily,  SEXP sZeroTol);

#endif
