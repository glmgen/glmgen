#include "utils.h"
#include "cs.h"

/* Calculates A*b + I (identity), for a square matrix
   A and scalar b */
cs * scalar_plus_eye (const cs * A, double b)
{
  int i;
  int j;
  cs * B;

  B = cs_spalloc(A->m, A->n, A->nzmax, 1, 0);

  for (j = 0; j < A->n; j++)
  {
    B->p[j] = A->p[j];
    for (i = A->p[j] ; i < A->p[j+1] ; i++)
    {
      if(A->i[i] == j)
      {
        B->x[i] = b * A->x[i] + 1;
      } else {
        B->x[i] = b * A->x[i];
      }
      B->i[i] = A->i[i];
    }
  }
  B->p[j] = A->p[j];

  return B;
}
