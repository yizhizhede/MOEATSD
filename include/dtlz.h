#ifndef _DTLZ_H 
#define _DTLZ_H

#include "problem.h"

Problem_t *DTLZ_new (char *title, int numObj, int numVar);
Matrix_t  *DTLZ_sample (int No, int numObj);

#endif
