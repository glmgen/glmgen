#include "cs.h"
#include <R_ext/Print.h>

/* print a sparse matrix; use %g for integers to avoid differences with csi */
csi cs_print (const cs *A, csi brief)
{
    csi p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { Rprintf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    if (nz < 0)
    {
        Rprintf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n]), cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            Rprintf ("    col %g : locations %g to %g\n", (double) j,
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                Rprintf ("      %g : %g\n", (double) (Ai [p]), Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { Rprintf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        Rprintf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            Rprintf ("    %g %g : %g\n", (double) (Ai [p]), (double) (Ap [p]),
                Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { Rprintf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}
