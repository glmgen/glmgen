#include "utils.h"
#include "cs.h"

/* Takes gqr structure from matrix A and solves Ax=b, overwriting the input b */
csi glmgen_qrsol (gqr * B, double * b)
{
  csi k;

  if(B->m < B->n) return(1); /* we only deal with m >= n case here */
  cs_ipvec (B->S->pinv, b, B->W, B->m) ;   /* x(0:m-1) = b(p(0:m-1) */
  for (k = 0 ; k < B->n ; k++)       /* apply Householder refl. to x */
  {
    cs_happly (B->N->L, k, B->N->B [k], B->W) ;
  }
  cs_usolve (B->N->U, B->W) ;           /* x = R\x */
  cs_ipvec (B->S->q, B->W, b, B->n) ;      /* b(q(0:n-1)) = x(0:n-1) */

  return (1) ;
}
