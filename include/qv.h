#ifndef _QV_H 
#define _QV_H

#include "problem.h"

Problem_t *QV_new (char *title, int numObj, int numVar);
Matrix_t  *QV_sample (int No, int numObj);

#endif
