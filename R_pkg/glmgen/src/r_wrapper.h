#ifndef RWRAPPER_H
#define RWRAPPER_H

double get_control_value(SEXP sControlList, const char * param_name, double param_default);

SEXP tf_R ( SEXP sY, SEXP sX, SEXP sN, SEXP sK, SEXP sFamily, SEXP sMethod, SEXP sMaxIter,
            SEXP sLamFlag, SEXP sObjFlag, SEXP sLambda, SEXP sNlambda, SEXP sLambdaMinRatio,
            SEXP sControl );

#endif