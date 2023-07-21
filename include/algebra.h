#ifndef _ALGEBRA_H
#define _ALGEBRA_H

#include "matrix.h"

int    	sum  (int *A, int n);
int    	isLinear (double *X, double *Y, int n);
double 	mean (double *A, int n);
double 	VAR  (double *A, int n);
double 	STD  (double *A, int n);
double 	norm (double *A, int n);
double 	sum  (double *A, int n);
double 	inner_product (double *A, double *B, int n);
double*	scale_product (double s, double *A, int n);
double*	vector_substract (double *A, double *B, int n);
double 	distance_p2l (double *p, double *line, int n);
double 	distance_p2p (double *p1, double *p2, int n);
double 	vector_angle (double *v1, double *v2, int n);
double 	distance_p2m (double *point, Matrix_t *M);
double 	distance_m2m (Matrix_t *M1, Matrix_t *M2);
double 	Lp_norm (double *A, int n, double p);

Matrix_t* 	scale_product (double s, Matrix_t *M);

void 	make_line (Matrix_t *X, double *Y, double *belta, double &delta);
void 	mean (double *X, int rowDim, int colDim, double *X_bar); 	// averay by colDim 

#endif
