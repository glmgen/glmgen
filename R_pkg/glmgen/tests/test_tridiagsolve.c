

#include "stdio.h"
#include "tridiagsolve.c"

int main() {
  int i;
	int  n = 4;
	double a[3] = { -1, -1, -1 };
	double b[4] = { 4,  4,  4,  4 };
	double c[3] = {-1, -1, -1 };
	double d[4] = { 5,  5, 10, 23 };
	double scratch[4] = { 0, 0, 0, 0 };
	// results    { 2,  3,  5, 7  }
	tridiagsolve(n,a,b,c,d,scratch);
	for (i = 0; i < n; i++) {
		printf( "%.0f\n", d[i] );
	}

	return 0;
}
