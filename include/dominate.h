#ifndef _DOMINATE_H
#define _DOMINATE_H

#include "matrix.h"
#include "link.h"

int isDominate   (double *a, double *b, int n);
int *isDominated (Matrix_t *A, Matrix_t *B);

int isDominate  (Vector_t *v1, Vector_t *v2);
int isDominated (Vector_t *v, Matrix_t *M);

List_t *ndSort (Matrix_t *M);
List_t *ndSort (Matrix_t *M, int minNum);


#endif
