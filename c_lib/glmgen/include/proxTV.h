
#ifndef PROXTV_H
#define PROXTV_H

/* Comparison tolerance */
#define EPSILON 1e-10
#define IS_ZERO(x) (x < EPSILON & x > -EPSILON)
#define IS_POSITIVE(x) (x > EPSILON)
#define IS_NEGATIVE(x) (x < -EPSILON)

int tautString_TV1(double *y,double lambda,double *x,int n);
int tautString_TV1_Weighted(double *y,double *lambda,double *x,int n);

#endif
