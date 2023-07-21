#include "dominate.h"
#include <stdlib.h>
#include <float.h>

int isDominate (double *a, double *b, int n) {
	int 	i;
	int 	positive = 0;
	int 	equals = 0;	
	double 	d;

	for (i=0; i<n; i++) {
		d = a[i] - b[i];
		if (d < 0)
			positive++;
		else if (0 == d)
			equals++;
		else
			return -1;
	}
	if (equals < n)
		return 1;
	else 
		return 0;
}

int *isDominated (Matrix_t *A, Matrix_t *B)
{
	int *flag = (int *)malloc (A->rowDim*sizeof (int));
	int i, j, col = A->colDim;

	for (i=0; i < A->rowDim; i++) {
		flag[i] = 0;
		for (j=0; j < B->rowDim; j++) {
			if (isDominate(B->elements+j*col, A->elements+i*col, col) == 1) {
				flag[i] = 1;
				break;
			}
		}
	}

	return flag;
}

int isDominate (Vector_t *v1, Vector_t *v2) {
	return isDominate (v1->elements, v2->elements, v1->dim);
}

int isDominated (Vector_t *v, Matrix_t *M) {

	int i, col = M->colDim;

	for (i=0; i < M->rowDim; i++) {
		if (isDominate(M->elements+i*col, v->elements, col) == 1) {
			return 1; 
		}
	}

	return 0;
}
