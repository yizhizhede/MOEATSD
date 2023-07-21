#ifndef _LANDSCAPE_H
#define _LANDSCAPE_H

// the landscape function of LSMOP
double Sphere (double *x, int n);
double Schwefel (double *x, int n);
double Rosenbrock (double *x, int n);
double Rastrigin (double *x, int n);
double Griewank (double *x, int n);
double Ackley (double *x, int n);

// the landscape function of DTLZ
double dtlz_g1 (double *x, int n);
double dtlz_g2 (double *x, int n);
double dtlz_g6 (double *x, int n);
double dtlz_g7 (double *x, int n);

#endif
