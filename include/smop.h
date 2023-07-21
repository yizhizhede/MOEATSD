#ifndef _SMOP_H 
#define _SMOP_H

#include "problem.h"

Problem_t *SMOP_new (char *title, int numObj, int numVar);
Matrix_t  *SMOP_sample (int No, int numObj);

#endif
