#include "stdio.h"
#include "tridiagsolve.c"

int main() {
  int i;
  int  n = 4;
  //double a[3] = { -1, -1, -1 };
  //double b[4] = { 4,  4,  4,  4 };
  //double c[3] = {-1, -1, -1 };
  //double d[4] = { 5,  5, 10, 23 };
  //double scratch[4] = { 0, 0, 0, 0 };
  double * a = (double *)malloc(3*sizeof(double)); a[0] = a[1] = a[2] = -1;
  double * b = (double *)malloc(4*sizeof(double)); b[0] = b[1] = b[2] = b[3] = 4;
  double * c = (double *)malloc(3*sizeof(double)); c[0] = c[1] = c[2] = -1;
  double * d = (double *)malloc(4*sizeof(double)); 
  d[0] = 5; d[1] = 5; d[2] = 10; d[3] = 23;

  double * scratch = (double *)malloc(4*sizeof(double)); 	

  // results    { 2,  3,  5, 7  }
  tridiagsolve(n,a,b,c,d,scratch);
  for (i = 0; i < n; i++) {
    printf( "%.0f\n", d[i] );
  }

  free(a); free(b); free(c); free(d); free(scratch);
  return 0;
}

// gcc -g -I../src -I../inst/include test_tridiagsolve.c
