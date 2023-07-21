#ifndef _NSGA3_H
#define _NSGA3_H

#include "population.h"
#include "problem.h"

Population_t* nsga3 (Problem_t *problem);
Population_t* nsga3_eliminate (Population_t *Q, Matrix_t *Z);

#endif
