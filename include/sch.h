#ifndef _SCH_H 
#define _SCH_H

#include "problem.h"

Problem_t *SCH_new (char *title, int numObj, int numVar);
Matrix_t  *SCH_sample (int No, int numObj);

#endif
