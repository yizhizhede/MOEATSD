#include "crowding.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double *crowding_distance (Matrix_t *M) {
	int 	i, j;
	int*	index = NULL;
	double*	I = (double *)malloc (M->rowDim*sizeof (double));
	double 	d = 0;

	memset (I, 0, M->rowDim*sizeof (double));
	for (j=0; j<M->colDim; j++) {
		index = sort (M, j);
		d = M->elements[index[M->rowDim-1]*M->colDim + j] - M->elements[index[0]*M->colDim +j];
		if (d < 1.0e-10) {
			free (index);
			continue;
		}

		I[index[0]] = I[index[M->rowDim-1]] = 1.0e+10;
		for (i=1; i<M->rowDim-1; i++) {
			I[index[i]] += (M->elements[index[i+1]*M->colDim+j] - M->elements[index[i-1]*M->colDim+j]) / d; 
		}
		free (index);
	}

	return I;
}

