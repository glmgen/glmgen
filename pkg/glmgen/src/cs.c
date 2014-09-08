// Code modified from: CSparse: a Concise Sparse Matrix package. VERSION 3.1.3
// Copyright (c) 2006-2014, Timothy A. Davis, Mar 26, 2014. Released under LGPL
// URL: http://www.suitesparse.com

#include "cs.h"

void *cs_malloc (csi n, size_t size)
{
    return (malloc (CS_MAX (n,1) * size)) ;
}

csn *cs_ndone (csn *N, cs *C, void *w, void *x, csi ok)
{
    cs_spfree (C) ;                     /* free temporary matrix */
    cs_free (w) ;                       /* free workspace */
    cs_free (x) ;
    return (ok ? N : cs_nfree (N)) ;    /* return result if OK, else free it */
}

cs *cs_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    cs_free (A->p) ;
    cs_free (A->i) ;
    cs_free (A->x) ;
    return ((cs *) cs_free (A)) ;   /* free the cs struct and return NULL */
}

void *cs_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet)
{
    cs *A = cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = cs_malloc (triplet ? nzmax : n+1, sizeof (csi)) ;
    A->i = cs_malloc (nzmax, sizeof (csi)) ;
    A->x = values ? cs_malloc (nzmax, sizeof (double)) : NULL ;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

void *cs_calloc (csi n, size_t size)
{
    return (calloc (CS_MAX (n,1), size)) ;
}

csi cs_happly (const cs *V, csi i, double beta, double *x)
{
    csi p, *Vp, *Vi ;
    double *Vx, tau = 0 ;
    if (!CS_CSC (V) || !x) return (0) ;     /* check inputs */
    Vp = V->p ; Vi = V->i ; Vx = V->x ;
    for (p = Vp [i] ; p < Vp [i+1] ; p++)   /* tau = v'*x */
    {
        tau += Vx [p] * x [Vi [p]] ;
    }
    tau *= beta ;                           /* tau = beta*(v'*x) */
    for (p = Vp [i] ; p < Vp [i+1] ; p++)   /* x = x - v*tau */
    {
        x [Vi [p]] -= Vx [p] * tau ;
    }
    return (1) ;
}

csi cs_scatter (const cs *A, csi j, double beta, csi *w, double *x, csi mark,
    cs *C, csi nz)
{
    csi i, p, *Ap, *Ai, *Ci ;
    double *Ax ;
    if (!CS_CSC (A) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}

double cs_house (double *x, double *beta, csi n)
{
    double s, sigma = 0 ;
    csi i ;
    if (!x || !beta) return (-1) ;          /* check inputs */
    for (i = 1 ; i < n ; i++) sigma += x [i] * x [i] ;
    if (sigma == 0)
    {
        s = fabs (x [0]) ;                  /* s = |x(0)| */
        (*beta) = (x [0] <= 0) ? 2 : 0 ;
        x [0] = 1 ;
    }
    else
    {
        s = sqrt (x [0] * x [0] + sigma) ;  /* s = norm (x) */
        x [0] = (x [0] <= 0) ? (x [0] - s) : (-sigma / (x [0] + s)) ;
        (*beta) = -1. / (s * x [0]) ;
    }
    return (s) ;
}

csn *cs_nfree (csn *N)
{
    if (!N) return (NULL) ;     /* do nothing if N already NULL */
    cs_spfree (N->L) ;
    cs_spfree (N->U) ;
    cs_free (N->pinv) ;
    cs_free (N->B) ;
    return ((csn *) cs_free (N)) ;  /* free the csn struct and return NULL */
}

csn *cs_qr (const cs *A, const css *S)
{
    double *Rx, *Vx, *Ax, *x,  *Beta ;
    csi i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, *s, *leftmost, *Ap, *Ai,
        *parent, *Rp, *Ri, *Vp, *Vi, *w, *pinv, *q ;
    cs *R, *V ;
    csn *N ;
    if (!CS_CSC (A) || !S) return (NULL) ;
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    q = S->q ; parent = S->parent ; pinv = S->pinv ; m2 = S->m2 ;
    vnz = S->lnz ; rnz = S->unz ; leftmost = S->leftmost ;
    w = cs_malloc (m2+n, sizeof (csi)) ;            /* get csi workspace */
    x = cs_malloc (m2, sizeof (double)) ;           /* get double workspace */
    N = cs_calloc (1, sizeof (csn)) ;               /* allocate result */
    if (!w || !x || !N) return (cs_ndone (N, NULL, w, x, 0)) ;
    s = w + m2 ;                                    /* s is size n */
    for (k = 0 ; k < m2 ; k++) x [k] = 0 ;          /* clear workspace x */
    N->L = V = cs_spalloc (m2, n, vnz, 1, 0) ;      /* allocate result V */
    N->U = R = cs_spalloc (m2, n, rnz, 1, 0) ;      /* allocate result R */
    N->B = Beta = cs_malloc (n, sizeof (double)) ;  /* allocate result Beta */
    if (!R || !V || !Beta) return (cs_ndone (N, NULL, w, x, 0)) ;
    Rp = R->p ; Ri = R->i ; Rx = R->x ;
    Vp = V->p ; Vi = V->i ; Vx = V->x ;
    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */
    rnz = 0 ; vnz = 0 ;
    for (k = 0 ; k < n ; k++)               /* compute V and R */
    {
        Rp [k] = rnz ;                      /* R(:,k) starts here */
        Vp [k] = p1 = vnz ;                 /* V(:,k) starts here */
        w [k] = k ;                         /* add V(k,k) to pattern of V */
        Vi [vnz++] = k ;
        top = n ;
        col = q ? q [k] : k ;
        for (p = Ap [col] ; p < Ap [col+1] ; p++)   /* find R(:,k) pattern */
        {
            i = leftmost [Ai [p]] ;         /* i = min(find(A(i,q))) */
            for (len = 0 ; w [i] != k ; i = parent [i]) /* traverse up to k */
            {
                s [len++] = i ;
                w [i] = k ;
            }
            while (len > 0) s [--top] = s [--len] ; /* push path on stack */
            i = pinv [Ai [p]] ;             /* i = permuted row of A(:,col) */
            x [i] = Ax [p] ;                /* x (i) = A(:,col) */
            if (i > k && w [i] < k)         /* pattern of V(:,k) = x (k+1:m) */
            {
                Vi [vnz++] = i ;            /* add i to pattern of V(:,k) */
                w [i] = k ;
            }
        }
        for (p = top ; p < n ; p++) /* for each i in pattern of R(:,k) */
        {
            i = s [p] ;                     /* R(i,k) is nonzero */
            cs_happly (V, i, Beta [i], x) ; /* apply (V(i),Beta(i)) to x */
            Ri [rnz] = i ;                  /* R(i,k) = x(i) */
            Rx [rnz++] = x [i] ;
            x [i] = 0 ;
            if (parent [i] == k) vnz = cs_scatter (V, i, 0, w, NULL, k, V, vnz);
        }
        for (p = p1 ; p < vnz ; p++)        /* gather V(:,k) = x */
        {
            Vx [p] = x [Vi [p]] ;
            x [Vi [p]] = 0 ;
        }
        Ri [rnz] = k ;                     /* R(k,k) = norm (x) */
        Rx [rnz++] = cs_house (Vx+p1, Beta+k, vnz-p1) ; /* [v,beta]=house(x) */
    }
    Rp [n] = rnz ;                          /* finalize R */
    Vp [n] = vnz ;                          /* finalize V */
    return (cs_ndone (N, NULL, w, x, 1)) ;  /* success */
}