#ifndef _UF_H 
#define _UF_H

#include "problem.h"

Problem_t *UF_new (char *title, int numObj, int numVar);
Matrix_t  *UF_sample (int No, int numObj);

#endif
