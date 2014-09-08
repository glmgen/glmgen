// Code modified from: CSparse: a Concise Sparse Matrix package. VERSION 3.1.3
// Copyright (c) 2006-2014, Timothy A. Davis, Mar 26, 2014. Released under LGPL
// URL: http://www.suitesparse.com

// These are the minimum typedef and functions from CSparse needed to execute
//   a sparse QR.

#ifndef CS_H
#define CS_H
#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#define CS_CSC(A) (A && (A->nz == -1))
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))

#ifndef csi
#define csi ptrdiff_t
#endif

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    csi nzmax ;     /* maximum number of entries */
    csi m ;         /* number of rows */
    csi n ;         /* number of columns */
    csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    csi nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs *L ;         /* L for LU and Cholesky, V for QR */
    cs *U ;         /* U for LU, R for QR, not used for Cholesky */
    csi *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    csi *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    csi *q ;        /* fill-reducing column permutation for LU and QR */
    csi *parent ;   /* elimination tree for Cholesky and QR */
    csi *cp ;       /* column pointers for Cholesky, row counts for QR */
    csi *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    csi m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
} css ;

void *cs_malloc (csi n, size_t size);
csn *cs_ndone (csn *N, cs *C, void *w, void *x, csi ok);
cs *cs_spfree (cs *A);
void *cs_free (void *p);
cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet);
void *cs_calloc (csi n, size_t size);
csi cs_happly (const cs *V, csi i, double beta, double *x);
csi cs_scatter (const cs *A, csi j, double beta, csi *w, double *x, csi mark,
    cs *C, csi nz);
double cs_house (double *x, double *beta, csi n);
csn *cs_nfree (csn *N);
csn *cs_qr (const cs *A, const css *S);

#endif
