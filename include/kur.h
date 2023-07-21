#ifndef _KUR_H 
#define _KUR_H

#include "problem.h"

Problem_t *KUR_new (char *title, int numObj, int numVar);
Matrix_t  *KUR_sample (int No, int numObj);

#endif
