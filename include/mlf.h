#ifndef _MLF_H 
#define _MLF_H

#include "problem.h"

Problem_t *MLF_new (char *title, int numObj, int numVar);
Matrix_t  *MLF_sample (int No, int numObj);

#endif
