#ifndef _WFG_H 
#define _WFG_H

#include "problem.h"
#include "matrix.h"

// genereate a wfg problem
Problem_t *WFG_new (char *title, int numObj, int numVar);

// a set of solution to WFG suite 
Matrix_t *WFG_sample (int No, int M);

#endif
