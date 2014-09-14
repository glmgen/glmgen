#ifndef UTILS_H
#define UTILS_H

double choose(int n, int k);
int max(int a, int b);
int min(int a, int b);
void soft_thresh(int n, double *y, double lam, double *beta);
gqr * glmgen_qr (const cs * A);
csi glmgen_qrsol (gqr * B, double * b);
csi glmgen_gqr_free (gqr * A);
cs * scalar_plus_eye (const cs * A, double b);

#endif


{
  gqr * B = malloc(1 * sizeof(cs));
  B->m = A->m;
  B->n = A->n;
  B->S = cs_sqr (1, A, 1) ;
  B->N = cs_qr (A, (B->S)) ;
  B->W = cs_calloc ((B->S) ? ((B->S))->m2 : 1, sizeof (double)) ;
  return (B);
}

/* Takes gqr structure from matrix A and solves Ax=b, overwriting the input b */
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

/* Free a constructed gqr struct */
{
  cs_sfree(A->S);
  cs_nfree(A->N);
  free(A->W);
  free(A);
  return(0);
}

/* Calculates A*b + I (identity), for a square matrix
   A and scalar b */
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