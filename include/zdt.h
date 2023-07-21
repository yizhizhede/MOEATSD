#ifndef _ZDT_H 
#define _ZDT_H

#include "problem.h"

Problem_t *ZDT_new (char *title, int numObj, int numVar);
Matrix_t  *ZDT_sample (int No, int numObj);

#endif
