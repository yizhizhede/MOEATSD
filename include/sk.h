#ifndef _SK_H 
#define _SK_H

#include "problem.h"

Problem_t *SK_new (char *title, int numObj, int numVar);
Matrix_t  *SK_sample (int No, int numObj);

#endif
