#ifndef _SNAPSHOT_H
#define _SNAPSHOT_H

#include "population.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

typedef struct snapshot_tag {
	char alg[32];
	Population_t *pop;	

	clock_t startTime;
	clock_t endTime;
     	double  runtime;
	
} snapshot_t;

void snapshot_init (char *alg, Population_t *pop);
void snapshot_click (Population_t *pop);

#endif
