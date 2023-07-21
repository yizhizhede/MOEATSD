#ifndef _ZLT_H 
#define _ZLT_H

#include "problem.h"

Problem_t *ZLT_new (char *title, int numObj, int numVar);
Matrix_t  *ZLT_sample (int No, int numObj);

#endif
