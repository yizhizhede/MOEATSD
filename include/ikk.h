#ifndef _IKK_H 
#define _IKK_H

#include "problem.h"

Problem_t *IKK_new (char *title, int numObj, int numVar);
Matrix_t  *IKK_sample (int No, int numObj);

#endif
