#ifndef _FA_H 
#define _FA_H

#include "problem.h"

Problem_t *FA_new (char *title, int numObj, int numVar);
Matrix_t  *FA_sample (int No, int numObj);

#endif
