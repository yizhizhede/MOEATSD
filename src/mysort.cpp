#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


static int maxDim = 0;
static int BYCOL = 0;
static int cmp_asc (const void *a, const void *b) {
	int 	i;
	double 	*aa = (double *)a;
	double 	*bb = (double *)b;
	double  t = 0;

	for (i=0; i<maxDim; i++) {
		t = aa[i] - bb[i];
		if (t < -DBL_EPSILON) return -1;
		if (t >  DBL_EPSILON) return 1;
	}
	return 0;
}

static int cmp_des (const void *a, const void *b) {
	int 	i;
	double 	*aa = (double *)a;
	double 	*bb = (double *)b;
	double  t=0;

	for (i=0; i<maxDim; i++) {
		t = aa[i] - bb[i];
		if ( t > DBL_EPSILON) return -1;
		if ( t < -DBL_EPSILON) return 1;
	}
	return 0;
}

static int cmp_asc_by (const void *a, const void *b) {
	int 	i = BYCOL;
	double 	*aa = (double *)a;
	double 	*bb = (double *)b;
	double  t;

	t = aa[i] - bb[i];
	if (t < -DBL_EPSILON) return -1;
	if (t > DBL_EPSILON) return 1;
	return 0;
}

static int cmp_des_by (const void *a, const void *b) {
	int 	i = BYCOL;
	double 	*aa = (double *)a;
	double 	*bb = (double *)b;
	double 	t;

	t = aa[i] - bb[i];
	if (t < -DBL_EPSILON) return 1;
	if (t > DBL_EPSILON) return -1;

	return 0;
}

int *sort (Matrix_t *M) {
	return sort (M, -1, "ASC");
}

int *sort (Matrix_t *M, int byCol) {
	return sort (M, byCol, "ASC");
}

int *sort (Matrix_t *M, const char* order) {
	return sort (M, -1, order);
}



int *sort (Matrix_t *M, int byCol, const char* order) {
	double*	buff = NULL;
	size_t 	size;
	int*	index = NULL;
	int	i;	

	maxDim = M->colDim + 1;
	BYCOL = byCol;

	size = M->rowDim * (M->colDim + 1) * sizeof (double);
	buff = (double *)malloc (size);
	
	size = M->colDim * sizeof (double);
	for (i=0; i<M->rowDim; i++) {
		memcpy (buff+i*(M->colDim + 1), M->elements+i*M->colDim, size);
		buff[i*(M->colDim + 1) + M->colDim] = 1.0 * i + 0.1;
	}
	
	size = (M->colDim + 1) * sizeof (double);
	if (byCol > -1) {
		if (strcmp (order,"DES") == 0 )
			qsort (buff, M->rowDim, size, cmp_des_by);
 		else 
			qsort (buff, M->rowDim, size, cmp_asc_by);
	} else {
		if (strcmp (order,"DES") == 0 )
			qsort (buff, M->rowDim, size, cmp_des);
 		else 
			qsort (buff, M->rowDim, size, cmp_asc);
	}


	size = M->rowDim * sizeof (int);
	index = (int *)malloc (size);
	for (i=0; i<M->rowDim; i++) {
		index[i] = (int)(buff[i*(M->colDim + 1) + M->colDim]);
	}
	
	free (buff);
	return index;
}
 

