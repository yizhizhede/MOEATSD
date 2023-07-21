#ifndef _SP_H 
#define _SP_H

#include "problem.h"

Problem_t *SP_new (char *title, int numObj, int numVar);
Matrix_t  *SP_sample (int No, int numObj);

#endif
