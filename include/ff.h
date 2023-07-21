#ifndef _FF_H 
#define _FF_H

#include "problem.h"

Problem_t *FF_new (char *title, int numObj, int numVar);
Matrix_t  *FF_sample (int No, int numObj);

#endif
