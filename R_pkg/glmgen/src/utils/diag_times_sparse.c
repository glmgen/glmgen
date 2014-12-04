#include "utils.h"
#include "cs.h"

/* Sets A = w*A */
void diag_times_sparse (const cs * A, double * w)
{
  int i;
  int j;

  for (j = 0; j < A->n; j++)
  {
    for (i = A->p[j] ; i < A->p[j+1] ; i++)
    {
      A->x[i] *= w[ A->i[i] ];
    }
  }

}
