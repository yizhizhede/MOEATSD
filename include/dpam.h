#ifndef _DPAM_H 
#define _DPAM_H

#include "problem.h"

Problem_t *DPAM_new (char *title, int numObj, int numVar);
Matrix_t  *DPAM_sample (int No, int numObj);

#endif
