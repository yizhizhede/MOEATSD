#ifndef _SSFYY_H 
#define _SSFYY_H

#include "problem.h"

Problem_t *SSFYY_new (char *title, int numObj, int numVar);
Matrix_t  *SSFYY_sample (int No, int numObj);

#endif
