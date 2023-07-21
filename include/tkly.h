#ifndef _TKLY_H 
#define _TKLY_H

#include "problem.h"

Problem_t *TKLY_new (char *title, int numObj, int numVar);
Matrix_t  *TKLY_sample (int No, int numObj);

#endif
