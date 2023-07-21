#ifndef _LMF_H 
#define _LMF_H

#include "problem.h"

Problem_t *LMF_new (char *title, int numObj, int numVar);
Matrix_t  *LMF_sample (int No, int numObj);

#endif
