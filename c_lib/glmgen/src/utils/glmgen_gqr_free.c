#include "utils.h"
#include "cs.h"

/* Free a constructed gqr struct */
csi glmgen_gqr_free (gqr * A)
{
  cs_sfree(A->S);
  cs_nfree(A->N);
  free(A->W);
  free(A);
  return(0);
}
