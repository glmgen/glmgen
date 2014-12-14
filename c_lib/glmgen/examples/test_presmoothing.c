#include "tf.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void test_presmoothing()
{
  double *x, *y, *w;
  double *xt, *yt, *wt;
  int n, nt;
  int k;
  int i;
  double delta;
  double x_cond = 1e15;

  srand(time(NULL));
  n = 10000; k = 2;

  x = (double*)malloc( n * sizeof(double) );
  y = (double*)malloc( n * sizeof(double) );
  w = (double*)malloc( n * sizeof(double) );

  x[0] = 0;
  for(i=1; i<n; i++) x[i] = x[i-1]+ (rand() % 100) * .02/(double)n;
  for(i=0; i<n; i++) y[i] = x[i] * x[i];
  for(i=0; i<n; i++) w[i] = 1;

  delta = 6 * (x[n-1] - x[0]) / n;
  
 
  xt = yt = wt = NULL;

  nt = 0;
  presmooth(x,y,w,n,k,&xt,&yt,&wt,&nt,x_cond);
  printf("delta = %g\n", delta); 
  printf("nt = %d\n", nt);
  //printf("x: "); for(i=0; i < nt; i++) printf("%g ", xt[i]); printf("\n\n");
  //printf("y: "); for(i=0; i < nt; i++) printf("%g ", yt[i]); printf("\n\n");
  //printf("w: "); for(i=0; i < nt; i++) printf("%g ", wt[i]); printf("\n\n");
 
  //printf("dx: ");
  double mindx = x[n-1] - x[0];
  for(i=0; i < nt -1; i++) {
    mindx = MIN(xt[i+1] - xt[i], mindx);
    if(mindx <= 0) printf("mindx <= 0 at i=%d\n", i);
    //printf("%g ", xt[i+1] - xt[i]);
  }
  printf("minimum dx = %g\n", mindx);
  printf("new range= %g, old range = %g\n", xt[nt-1]-xt[0], x[n-1]-x[0]);
  
  free(x); free(y); free(w);
  free(xt); free(yt); free(wt);

}
int main()
{
  test_presmoothing();
  return 0;
}

