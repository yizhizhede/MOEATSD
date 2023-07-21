#ifndef _SHAPE_H
#define _SHAPE_H

#include "matrix.h"

void H1 (double *x, double *g, double *f, int M);
void H2 (double *x, double *g, double *f, int M);
void H3 (double *x, double *g, double *f, int M);

Matrix_t *H1_sample (int M, int H);
Matrix_t *H2_sample (int M, int H);
Matrix_t *H3_sample (int M, int H);
Matrix_t *H4_sample (int M, int H);

Matrix_t *Grid_sample (int Dim, int H);

void linear (double *x, double *f, int M);
void convex (double *x, double *f, int M);
void concave (double *x, double *f, int M);
void mixedM (double *x, double *f, int M);
void discM (double *x, double *f, int M);

int conb (int n, int m);	// conbination: to select m individuals from n $C_n^m$
int grid (int dim, int H);	

#endif
