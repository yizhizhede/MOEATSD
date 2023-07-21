#include <stdio.h>
#include <stdlib.h>

#include "hv.h"
#include "matrix.h"

/*
static double intersection (Matrix_t *M) {
	double 	I = 1.0;
	int 	i = 0;
	Matrix_t *p = NULL;

	p = Matrix_min (M);	

	for (i=0; i < p->colDim; i++)
		I *= p->elements[i];
	Matrix_free (&p);
	return I;
}
*/

static double intersectionII (Matrix_t *M) {
	double 	I = 1.0;
	int 	i = 0;
	Matrix_t *p = NULL;

	p = Matrix_max (M);	

	for (i=0; i < p->colDim; i++)
		I *=(1.1 - p->elements[i]);
	Matrix_free (&p);
	return I;
}


double hv_iea (Matrix_t *M)
{
	unsigned int NO;
	int queue[100], tail = 0;
	int i;
	double hv = 0, I;
	Matrix_t *subM = NULL;

	if (M->rowDim >30) {
		return -1.0;
	}
	
	for (NO = ((1<<M->rowDim) - 1); NO > 0; NO--) {
		tail = 0;
		for (i=0; i<M->rowDim; i++) {
			if (((NO >> i) & 1) == 1) {
				queue[tail++] = i;
			}
		}

		subM = Matrix_sub (M, queue, tail);
		// I = intersection (subM);	
		I = intersectionII (subM);	
		Matrix_free (&subM);

		if (tail % 2 == 0) {
			hv -= I;
		} else {
			hv += I;
		}	
	}
	return hv;
}

