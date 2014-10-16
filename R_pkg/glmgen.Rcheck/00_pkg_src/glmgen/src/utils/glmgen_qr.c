#include "utils.h"
#include "cs.h"

/* Calculates symbolic and numeric qr for a given matrix A */
gqr * glmgen_qr (const cs * A)
{
  gqr * B = malloc(1 * sizeof(cs));
  B->m = A->m;
  B->n = A->n;
  B->S = cs_sqr (1, A, 1) ;
  B->N = cs_qr (A, (B->S)) ;
  B->W = cs_calloc ((B->S) ? ((B->S))->m2 : 1, sizeof (double)) ;
  return (B);
}
