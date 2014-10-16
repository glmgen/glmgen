#include "utils.h"
#include "cs.h"

/* Takes gqr structure from matrix A and solves Ax=b, overwriting the input b */
csi glmgen_qrsol (gqr * A, double * b)
{
  csi k;

  if(A->m < A->n) return(1); /* we only deal with m >= n case here */
  cs_ipvec (A->S->pinv, b, A->W, A->m) ;   /* x(0:m-1) = b(p(0:m-1) */
  for (k = 0 ; k < A->n ; k++)       /* apply Householder refl. to x */
  {
    cs_happly (A->N->L, k, A->N->B [k], A->W) ;
  }
  cs_usolve (A->N->U, A->W) ;           /* x = R\x */
  cs_ipvec (A->S->q, A->W, b, A->n) ;      /* b(q(0:n-1)) = x(0:n-1) */

  return (1) ;
}
