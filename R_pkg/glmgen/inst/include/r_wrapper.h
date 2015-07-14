#ifndef RWRAPPER_H
#define RWRAPPER_H

#include "utils.h"
#include "cs.h"

double get_control_value(SEXP sControlList, const char * param_name);
cs * dgTMatrix_to_cs (SEXP input);

SEXP thin_R (SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sControl);
SEXP tf_R (SEXP sY, SEXP sX, SEXP sW, SEXP sN, SEXP sK, SEXP sFamily,
    SEXP sMethod, SEXP sLamFlag, SEXP sLambda, SEXP sNlambda,
    SEXP sLambdaMinRatio, SEXP sControl );
SEXP tf_predict_R (SEXP sBeta, SEXP sX, SEXP sN, SEXP sK, SEXP sX0, SEXP sN0,
    SEXP sNLambda, SEXP sFamily,  SEXP sZeroTol);
SEXP lattice_R (SEXP sY, SEXP sW, SEXP sLambda, SEXP sRho, SEXP sEps,
                SEXP sMaxiter, SEXP sVerbose, SEXP sNaflag,
                SEXP sLatticeType, SEXP sE, SEXP sC, SEXP sBeta0);

#endif
