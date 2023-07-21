#ifndef _VFM_H 
#define _VFM_H

#include "problem.h"

Problem_t *VFM_new (char *title, int numObj, int numVar);
Matrix_t  *VFM_sample (int No, int numObj);

#endif
