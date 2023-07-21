#ifndef _LSMOP_H 
#define _LSMOP_H

#include "problem.h"

Problem_t *LSMOP_new (char *title, int numObj, int numVar);
Matrix_t  *LSMOP_sample (int No, int numObj);

#endif
