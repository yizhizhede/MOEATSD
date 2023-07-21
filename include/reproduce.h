#ifndef _REPRODUCE_H
#define _REPRODUCE_H

#include "population.h"

// randomly generate individual
Population_t *Reproduce_one_by_random (Population_t *P);	// create one individual
Population_t *Reproduce_equ_by_random (Population_t *P);	// create n (PopSize)  individuals

// generate individual by tournament
// That means the population is sorted, It shoud be used together with a rank routine
Population_t *Reproduce_equ_by_tournament (Population_t *P);	

#endif
