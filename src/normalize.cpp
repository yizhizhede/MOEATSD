#include "normalize.h"
#include "algebra.h"

Matrix_t *normalize (Matrix_t *M) {
	return Matrix_norm (M);
}

Matrix_t *normalize_sphere  (Matrix_t *M) {
	Matrix_t *dupM = Matrix_dup (M);
	int i, j, k;
	double r;

	for (i=0, k=0; i<M->rowDim; i++) {
		r = norm (M->elements+i*M->colDim, M->colDim);
		for (j=0; j<M->colDim; j++, k++)
			dupM->elements[k] /= (r > 0 ? r : 1 );
	}
	return dupM;	
}

