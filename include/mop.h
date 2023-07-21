#ifndef _MOP_H
#define _MOP_H

#include "problem.h"
#include "matrix.h"

// genereate a problem
Problem_t *MOP_new (char *title, int numObj, int numVar);

// a set of solution to WFG suite 
Matrix_t *MOP_sample (int No, int M);

#endif
