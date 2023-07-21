#ifndef _JOS_H 
#define _JOS_H

#include "problem.h"

Problem_t *JOS_new (char *title, int numObj, int numVar);
Matrix_t  *JOS_sample (int No, int numObj);

#endif
