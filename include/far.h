#ifndef _FAR_H 
#define _FAR_H

#include "problem.h"

Problem_t *FAR_new (char *title, int numObj, int numVar);
Matrix_t  *FAR_sample (int No, int numObj);

#endif
