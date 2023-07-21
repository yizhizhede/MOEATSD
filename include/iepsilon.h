#ifndef _IEPSILON_H
#define _IEPSILON_H

#include "matrix.h"

double Iepsilon (double *y1, double *y2, int m);
double Iepsilon (Matrix_t *A, Matrix_t *B);

double Iepsilon_plus (double *y1, double *y2, int m);
double Iepsilon_plus (Matrix_t *A, Matrix_t *B);

double *Iepsilon_fitness (Matrix_t *A);
void    Ipesilon_update  (Matrix_t *A, double *fitness, int *valid, int r);


#endif
