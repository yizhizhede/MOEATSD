#ifndef _CONTROLLER_H
#define _CONTROLLER_H

#include "problem.h"
#include "individual.h"

#define CONTROL_RETURN 			0
#define CONTROL_REPRODUCE		1
#define CONTROL_ELIMINATE		2	
#define CONTROL_REPRODUCE_ELIMINATE	3

// 
void control_return 	(Population20_t *pop);
void control_reproduce 	(Population20_t *pop);
void control_eliminate 	(Population20_t *pop);

//
int Controller (Problem_t *problem, Population20_t *pop);

#endif
