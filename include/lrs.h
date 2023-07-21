#ifndef _LRS_H 
#define _LRS_H

#include "problem.h"

Problem_t *LRS_new (char *title, int numObj, int numVar);
Matrix_t  *LRS_sample (int No, int numObj);

#endif
