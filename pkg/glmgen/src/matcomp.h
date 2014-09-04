#ifndef MATCOMP_H
#define MATCOMP_H

void d(int n, int k, double *drow, double *a, double *x);
void dt(int n, int k, double *drow, double *a, double *x);
void builddrow(int k, double *drow);
void builddmat(int n, int k, double *x, double *dmat);
void builddtd(int n, int k, double *drow, double *dtd);
void buildddt(int n, int k, double *drow, double *ddt);
void h(int n, int k, double *a, double *x);
void ht(int n, int k, double *a, double *x);
void l2d(int n, int k, double *a, double *x);
void solvetd1(int n, double r, double *d, double *bp, double *dp, double *x);
void solvetd2(int n, double *a, double *b, double *c, double *d, double *bp,
        double *dp, double *x);
void buildl(int n, int k, double *l);
void buildltl(int n, int k, double *l, double *ltl);
void solveddt(int n, int k, double *b, double *a, double *c, double *l,
        double *x);

#endif