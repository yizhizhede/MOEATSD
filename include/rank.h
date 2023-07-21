#ifndef _RANK_H
#define _RANK_H

#include "population.h"
#include "matrix.h"


int *rank_by_crowding (Matrix_t *M);		// non-dominated sort + crowding distance
int *rank_by_crowding (Matrix_t *M, int T);	// non-dominated sort + crowding distance

Population_t *rank_by_crowding (Population_t *pop);
Population_t *rank_by_crowding (Population_t *pop, int T);

int *rank_by_density (Matrix_t *M); 		// non-dominated sort + Shift-based density estimation (SDE)

#endif
