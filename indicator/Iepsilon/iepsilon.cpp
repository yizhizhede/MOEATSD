#include "iepsilon.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

double Iepsilon (double *y1, double *y2, int m) {
	int i;
	double xValue=0;
	double t;

	for (i=0; i<m; i++) {
		t = (y2[i] != 0)? y1[i]/y2[i]:y1[i]*1.0e+10;
		if (t > xValue) {
			xValue = t;
		}
	}
	return xValue;
}

double Iepsilon (Matrix_t *A, Matrix_t *B) {
	double *z1 = A->elements;
	double *z2 = B->elements;
	int i, j, k;
	double xValue1,	xValue2, xValue3, t;

	// max
	for (i=0, xValue1=-1.0e+10; i<B->rowDim; i++, z2+=B->colDim) {
		// min
		for (j=0, xValue2=1.0e+10; j<A->rowDim; j++, z1+=A->colDim) {
			// max 
			for (k=0, xValue3 = -1.0e+10; k<A->colDim; k++) {
				t = (z2[i] != 0)? z1[i]/z2[i]:z1[i]*1.0e+10;
				if (t > xValue3) {
					xValue3 = t;
				}
			}
			if (xValue3 < xValue2) {
				xValue2 = xValue3;
			}
		}
		if (xValue2 > xValue1) {
			xValue1 = xValue2;
		}
	}
	return xValue1;
}

double Iepsilon_plus (double *y1, double *y2, int m) {
	int i;
	double xValue=y1[0] -y2[0];
	double t;

	for (i=1; i<m; i++) {
		t =  y1[i]-y2[i];
		if (t > xValue) {
			xValue = t;
		}
	}
	return xValue;
}

double Iepsilon_plus (Matrix_t *A, Matrix_t *B) {
	double *z1 = A->elements;
	double *z2 = B->elements;
	int i, j, k;
	double xValue1,	xValue2, xValue3, t;

	// max
	for (i=0, xValue1=-1.0e+10; i<B->rowDim; i++, z2+=B->colDim) {
		// min
		for (j=0, xValue2=1.0e+10; j<A->rowDim; j++, z1+=A->colDim) {
			// max 
			for (k=0, xValue3 = -1.0e+10; k<A->colDim; k++) {
				t = z1[k] - z2[k];
				if (t > xValue3) {
					xValue3 = t;
				}
			}
			if (xValue3 < xValue2) {
				xValue2 = xValue3;
			}
		}
		if (xValue2 > xValue1) {
			xValue1 = xValue2;
		}
	}
	return xValue1;
}

double *Iepsilon_fitness (Matrix_t *A) {
	Matrix_t *normA = Matrix_norm (A);
	int i, j;
	size_t size = A->rowDim * sizeof (double);
	double *fitness = (double *)malloc (size);
	double *x1 = NULL, *x2=NULL;

	for (i=0; i<A->rowDim; i++) {
		x1 = normA->elements + i*A->colDim;
		for (j=0, fitness[i]=0.0; j<A->rowDim; j++) {
			if (i==j) continue;
			x2 = normA->elements + j*A->colDim;
			fitness[i] += -exp(-Iepsilon_plus (x2, x1, A->colDim)/0.05);
		}
	}

	Matrix_free (&normA);
	return fitness;
}

void Ipesilon_update  (Matrix_t *A, double *fitness, int *valid, int r) {
	Matrix_t *normA = Matrix_norm (A);
	int i;
	double *x1 = NULL, *x2=NULL;

	for (i=0; i<A->rowDim; i++) if (valid[i]) {
		x1 = normA->elements + i*A->colDim;
		x2 = normA->elements + r*A->colDim;
		fitness[i] += exp(-Iepsilon_plus (x2, x1, A->colDim)/0.05);
	}

	Matrix_free (&normA);
}
