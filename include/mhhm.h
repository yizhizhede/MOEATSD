#ifndef _MHHM_H 
#define _MHHM_H

#include "problem.h"

Problem_t *MHHM_new (char *title, int numObj, int numVar);
Matrix_t  *MHHM_sample (int No, int numObj);

#endif
