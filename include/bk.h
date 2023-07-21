#ifndef _BK_H 
#define _BK_H

#include "problem.h"

Problem_t *BK_new (char *title, int numObj, int numVar);
Matrix_t  *BK_sample (int No, int numObj);

#endif
