#include <math.h>
#include <stdlib.h>
#include <R_ext/Lapack.h>

#include "matcomp.h"
#include "utils.h"

// Matrix multiplication by Dk, the kth order difference
// operator, i.e., compute x = Dk*a
// Parameters:
// - n is the number of columns in Dk
// - k is the difference order, so Dk has n-k rows
// - drow is the row pattern of Dk, precomputed
// - a is the vector being multiplied, length n
// - x is where we put the answer, length n-k

void d(int n, int k, double *drow, double *a, double *x)
{
  int i;
  int j;

  for (i=0; i<n-k; i++)
  {
    x[i] = 0;
    for (j=0; j<k+1; j++)
    {
      x[i] += drow[j]*a[i+j];
    }
  }
}

// Matrix multiplication by Dk^T, the transpose of the kth
// order difference operator, i.e., compute x = Dk^T*a
// Parameters:
// - n is the number of columns in Dk
// - k is the difference order, so Dk has n-k rows
// - drow is the row pattern of Dk, precomputed
// - a is the vector being multiplied, length n-k
// - x is where we put the answer, length n

void dt(int n, int k, double *drow, double *a, double *x)
{
  int i;
  int j;

  for (i=0; i<n; i++)
  {
    x[i] = 0;
    for (j=min(i,k); j>=max(i-n+k+1,0); j--)
    {
      x[i] += drow[j]*a[i-j];
    }
  }
}

// Build a row (with zeros removed) of the matrix Dk, the
// kth order difference operator
void builddrow(int k, double *drow)
{
  int i;

  for (i=0; i<k+1; i++)
  {
    drow[i] = (1-2*((i+k+1)%2==0)) * choose(k,i);
  }
}

//INCOMPLETE
void builddmat(int n, int k, double *x, double *dmat)
{
  int i;
  int j;
  int row;
  int col;
  double * tmp;

  tmp = (double*)malloc(sizeof(double)*(n*(n-1)));

  for (i=0; i<n-1; i++)
  {
      tmp[i*n+i] = -1;
      tmp[i*n+i+1] = 1;
  }
  //tmp should now be D

  for (j=1; j<=k; j++) //perform k-1 matrix-multiplicaions
  {
    for (row=0; row<n-j-1; row++) //form the matrix row by row
    {
      for(col=row; col<row+j+2; col++) //for each row, update only few cols
      {
        //tmp[row*(n-j)+col] = -1.0 * tmp[row*(n-j) + col] / (x[row+j]-x[row]) + 1.0 * [(row+1)*(n-j) + col] / (x[row+j+1] - x[row+1]);
        //multiplying D by D^j
      }
    }
    //tmp should now be D^{j+1}
  }
  //tmp should now be D^{k+1}

  for(i=0;i< n;i++)
  {
    for(j=0; j< n-k; j++)
    {
      dmat[i*n+j] = tmp[i*n+j];
    }
  }

  free(tmp);
}

// Build the matrix Dk*Dk^T, where Dk is the kth order
// difference operator
// Stored in "L" packed representation for symmetric band
// matrix (in *column-major* order! for Fortran)
// Parameters:
// - n is the number of columns in Dk
// - k is the difference order, so Dk has n-k rows
// - drow is the row pattern of Dk, precomputed
// - dtd is where we put the answer, containing only the
//   lower triangular part, and in column-major order

void buildddt(int n, int k, double *drow, double *ddt)
{
  // Build one column of ddt (since we're using packed
  // representation, this is just diagonal and down)
  int i;
  int j;
  int r;
  double * tmp;

  tmp = (double*) malloc(sizeof(double)*(k+1));
  for (i=0; i<k+1; i++)
  {
    tmp[i] = 0;
    for (r=0; r<k+1-i; r++)
    {
      tmp[i] += drow[i+r]*drow[r];
    }
  }
  // Populate ddt appropriately
  for (i=0; i<k+1; i++)
  {
    for (j=0; j<n-k; j++)
    {
      ddt[i+j*(k+1)] = tmp[i];
    }
  }

  free(tmp);
}

// Build the matrix Dk^T*Dk, where Dk is the kth order
// difference operator
// Stored in "L" packed representation for symmetric band
// matrix (in *column-major* order! for Fortran)
// Parameters:
// - n is the number of columns in Dk
// - k is the difference order, so Dk has n-k rows
// - drow is the row pattern of Dk, precomputed
// - dtd is where we put the answer, containing only the
//   lower triangular part, and in column-major order

void builddtd(int n, int k, double *drow, double *dtd)
{
  int i;
  int l;
  int r;
  double t;

  for (i=0; i<=k; i++)
  {
    for (l=i; l<n; l++)
    {
      t = 0;
      for (r=min(l,k); r>=max(l-n+k+1,i); r--)
      {
        t += drow[r]*drow[r-i];
      }
      dtd[i+(l-i)*(k+1)] = t;
    }
  }

  // // Populate first k and last k columns of dtd
  // // First k
  // for (int i=0; i<k+1; i++) {
  //   for (int l=i; l<i+k; l++) {
  //     dtd[i+(l-i)*(k+1)] = 0;
  //     for (int r=min(l,k); r>=max(l-n+k+1,i); r--) {
  //  dtd[i+(l-i)*(k+1)] += drow[r]*drow[r-i];
  //     }
  //   }
  //   // Last k
  //   for (int l=i+n-k; l<i+n; l++) {
  //     dtd[i+(l-i)*(k+1)] = 0;
  //     for (int r=min(l,k); r>=max(l-n+k+1,i); r--) {
  //  dtd[i+(l-i)*(k+1)] += drow[r]*drow[r-i];
  //     }
  //   }
  // }

  // // Now build middle n-2k columns. Since they're all the
  // // same, just build one such column, then populate dtd
  // double *tmp = (double*)malloc(sizeof(double)*(k+1));
  // for (int i=0; i<k+1; i++) {
  //   tmp[i] = 0;
  //   for (int r=0; r<k+1-i; r++) {
  //     tmp[i] += drow[i+r]*drow[r];
  //   }
  // }
  // for (int i=0; i<k+1; i++) {
  //   for (int j=k; j<n-k; j++) {
  //     dtd[i+j*(k+1)] = tmp[i];
  //   }
  // }
  // free(tmp);
}

// Matrix multiplication by Hk, the kth order falling
// factorial basis matrix, i.e., compute x = Hk*a
// Parameters:
// - n is the number of input points
// - k is the polynomial order
// - a is the vector being multiplied, length n
// - x is where we put the answer, length n

void h(int n, int k, double *a, double *x)
{
  int i;
  int j;

  // Initialize x to be a
  for (i=0; i<n; i++)
  {
    x[i] = a[i];
  }

  // Take k+1 growing cumulative sums
  for (i=0; i<k+1; i++)
  {
    for (j=k+1-i; j<n; j++)
    {
      x[j] += x[j-1];
    }
  }
}

// Matrix multiplication by Hk^T, the transpose of the
// kth order falling factorial basis matrix, i.e., compute
// x = Hk^T*a
// Parameters:
// - n is the number of input points
// - k is the polynomial order
// - a is the vector being multiplied, length n
// - x is where we put the answer, length n

void ht(int n, int k, double *a, double *x)
{
  int i;
  int j;

  // Initialize x to be a
  for (i=0; i<n; i++) {
    x[i] = a[i];
  }

  // Take k+1 shrinking cumulative sums
  for (i=0; i<k+1; i++)
  {
    for (j=n-2; j>=i; j--)
    {
      x[j] += x[j+1];
    }
  }
}

// Matrix multiplication by L2*D, where L2 denotes the last
// n-k-1 columns of Hk, and D denotes D(k+1), i.e., compute
// x = L2*D*a
// Parameters:
// - n is the number of input points
// - k is the polynomial order
// - a is the vector being multiplied, length n
// - x is where we put the answer, length n

void l2d(int n, int k, double *a, double *x)
{
  int i;
  double * drow;

  drow = (double*)malloc((k+2)*sizeof(double));

  builddrow(k+1,drow);
  d(n,k+1,drow,a,x+k+1); // Compute D*a

  // Compute L2*x
  for (i=0; i<k+1; i++) x[i] = 0;
  h(n,k,x,x);
  for (i=0; i<k+1; i++) x[i] = 0;

  free(drow);
}

// Tridiagonal matrix algorithm, for the special matrix:
//  1+r   -r
//   -r   1+2r  -r
//        -r    1+2r  -r
//                    ...
//                    -r   1+2r  -r
//                         -r    1+r
// Solves Tx = d. Parameters:
// - n is the size of the (square) tridiagonal matrix (we need n>=2)
// - r is the multiplier demonstrated above
// - d is the right hand side, length n
// - bp is a buffer, length n
// - dp is a buffer, length n
// - x is where we put the answer, length n

void solvetd1(int n, double r, double *d, double *bp, double *dp, double *x)
{
  int i;
  double m;

  bp[0] = 1+r;
  dp[0] = d[0];
  for (i=1; i<n-1; i++)
  {
    m = -r/bp[i-1];
    bp[i] = 1+2*r + m*r;
    dp[i] = d[i] - m*dp[i-1];
  }
  m = -r/bp[n-2];
  bp[n-1] = 1+r + m*r;
  dp[n-1] = d[n-1] - m*dp[n-2];

  x[n-1] = dp[n-1]/bp[n-1];
  for (i=n-2; i>=0; i--)
  {
    x[i] = (dp[i]+r*x[i+1])/bp[i];
  }
}

// Tridiagonal matrix algorithm
// See http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
// Solves Tx = d. Parameters:
// - n is the size of the (square) tridiagonal matrix (we need n>=1)
// - a is the subdiagonal, length n-1
// - b is the diagonal, length n
// - c is the superdiagonal, length n-1
// - d is the right-hand side, length n
// - bp is a buffer, length n
// - dp is a buffer, length n
// - x is where we put the solution, length n

void solvetd2(int n, double *a, double *b, double *c, double *d, double *bp,
        double *dp, double *x)
{
  int i;
  double m;

  bp[0] = b[0];
  dp[0] = d[0];
  for (i=1; i<n; i++)
  {
    m = a[i-1]/bp[i-1];
    bp[i] = b[i] - m*c[i-1];
    dp[i] = d[i] - m*dp[i-1];
  }

  x[n-1] = dp[n-1]/bp[n-1];
  for (i=n-2; i>=0; i--)
  {
    x[i] = (dp[i]-c[i]*x[i+1])/bp[i];
  }
}

void buildl(int n, int k, double *l)
{
  int i;
  int j;

  // The first column is all 1s
  for (i=0; i<n; i++) l[i] = 1;
  // The next k columns are given by cumulative sums
  for (j=1; j<k+1; j++)
  {
    for (i=0; i<n; i++)
    {
      if (i<j) l[i+j*n] = 0;
      else l[i+j*n] = l[i-1+j*n] + l[i-1+(j-1)*n];
    }
  }
}

void buildltl(int n, int k, double *l, double *ltl)
{
  int i;
  int j;
  int r;
  long double t;

  // Populate L
  // The first column is all 1s
  for (i=0; i<n; i++) l[i] = 1;
  // The next k columns are given by cumulative sums
  for (j=1; j<k+1; j++)
  {
    for (i=0; i<n; i++)
    {
      if (i<j) l[i+j*n] = 0;
      else l[i+j*n] = l[i-1+j*n] + l[i-1+(j-1)*n];
    }
  }

  // Compute L^T*L
  for (i=0; i<k+1; i++)
  {
    for (j=0; j<k+1; j++)
    {
      t = 0;
      for (r=0; r<n; r++)
      {
        t += l[r+i*n]*l[r+j*n];
      }
      ltl[i+j*(k+1)] = t;
    }
  }
}

void solveddt(int n, int k, double *b, double *a, double *c, double *l,
        double *x)
{
  int i;
  int nr;
  int nc=k+1;
  int nrhs=1;
  int nw=n*(k+1);
  int info=0;
  char trans;
  double * wrk;

  wrk = (double*) malloc(n*(k+1)*sizeof(double));

  // Populate a
  for (i=0; i<k+1; i++) a[i] = 0;
  for (i=k+1; i<n; i++) a[i] = b[i-k-1];

  // Compute Hk*a, zero out first k+1 components
  h(n,k,a,c);
  for (i=0; i<k+1; i++) c[i] = 0;

  // Compute (L^T*L)^{-1}*L^T*c, where L denotes
  // the first k+1 columns of Hk
  trans = 'n';
  nr=n;
  nc=k+1;
  nrhs=1;
  nw=n*(k+1);
  info=0;

  F77_CALL(dgels)(&trans,&nr,&nc,&nrhs,l,&nr,c,&nr,wrk,&nw,&info);

  // Concatenate a = (-c,b), and multiply Hk^T*Hk*a
  for (i=0; i<k+1; i++) a[i] = -c[i];
  for (i=k+1; i<n; i++) a[i] = b[i-k-1];
  h(n,k,a,c);
  ht(n,k,c,a);

  // Take the last n-k-1 components of a as the answer
  for (int i=0; i<n-k-1; i++) x[i] = a[i+k+1];

  free(wrk);

  // // Take the last n-k-1 components of c as the answer
  // for (int i=0; i<n-k-1; i++) x[i] = c[i+k+1];

  // // Compute Hk^T*c
  // ht(n,k,c,a);

  // // Compute (L^T*L)^{-1}*a, where L denotes
  // // the first k+1 columns of Hk
  // int m=k+1, nrhs=1, info=0;
  // F77_CALL(dpotrs)(&uplo,&m,&nrhs,ltl,&m,a,&m,&info);

  // // Concatenate c = (-a,b), and multiply Hk^T*Hk*c
  // for (int i=0; i<k+1; i++) c[i] = -a[i];
  // for (int i=k+1; i<n; i++) c[i] = b[i-k-1];
  // h(n,k,c,a);
  // ht(n,k,a,c);

  // // Take the last n-k-1 components of c as the answer
  // for (int i=0; i<n-k-1; i++) x[i] = c[i+k+1];
}
