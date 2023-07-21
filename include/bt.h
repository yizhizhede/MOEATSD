#ifndef _BT_H 
#define _BT_H

#include "problem.h"

Problem_t *BT_new (char *title, int numObj, int numVar);
Matrix_t  *BT_sample (int No, int numObj);

#endif
