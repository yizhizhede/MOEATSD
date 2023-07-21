#ifndef _IBEA_H 
#define _IBEA_H

#include "population.h"
#include "problem.h"

Population_t* ibea (Problem_t *problem);

Matrix_t* ibea_fitness_assign (Matrix_t *Obj);
void      ibea_fitness_update (Matrix_t *Obj, int r);

Population_t* ibea_environmental_select (Population_t* Q, int nPop);

#endif
