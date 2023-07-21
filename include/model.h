#ifndef _MODEL_H
#define _MODEL_H

#include "population.h"

void 	model_update (Population_t *pop);
double 	model_getValue (Population_t *pop, int depV, double* offset);

#endif
