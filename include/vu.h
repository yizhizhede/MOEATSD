#ifndef _VU_H 
#define _VU_H

#include "problem.h"

Problem_t *VU_new (char *title, int numObj, int numVar);
Matrix_t  *VU_sample (int No, int numObj);

#endif
