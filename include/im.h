#ifndef _IM_H 
#define _IM_H

#include "problem.h"

Problem_t *IM_new (char *title, int numObj, int numVar);
Matrix_t  *IM_sample (int No, int numObj);

#endif
