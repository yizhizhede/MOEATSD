#ifndef _NNTP_H 
#define _NNTP_H

#include "problem.h"

Problem_t *NNTP_new (char *title, int numObj, int numVar);
Matrix_t  *NNTP_sample (int No, int numObj);

#endif
