#ifndef _FES_H 
#define _FES_H

#include "problem.h"

Problem_t *FES_new (char *title, int numObj, int numVar);
Matrix_t  *FES_sample (int No, int numObj);

#endif
