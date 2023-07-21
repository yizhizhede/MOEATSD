#ifndef _LE_H 
#define _LE_H

#include "problem.h"

Problem_t *LE_new (char *title, int numObj, int numVar);
Matrix_t  *LE_sample (int No, int numObj);

#endif
