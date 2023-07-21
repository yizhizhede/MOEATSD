#ifndef _LOGS_H
#define _LOGS_H

#include "matrix.h"

void log_string (char *fn, char *str);
void log_matrix (char *fn, Matrix_t *M);
void log_vector (char *fn, int *buf, int len);
void log_vector (char *fn, double *buf, int len);

#endif
