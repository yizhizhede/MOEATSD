#ifndef _LTDZ_H 
#define _LTDZ_H

#include "problem.h"

Problem_t *LTDZ_new (char *title, int numObj, int numVar);
Matrix_t  *LTDZ_sample (int No, int numObj);

#endif
