#ifndef _DGO_H 
#define _DGO_H

#include "problem.h"

Problem_t *DGO_new (char *title, int numObj, int numVar);
Matrix_t  *DGO_sample (int No, int numObj);

#endif
