#include "linkage.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846264338327950288419716939937510

void L1 (double *x, int n, Problem_t *problem) {
	int 	i=0;
	int 	numObj = problem->numObj;
	int 	numVar = problem->numVar;
	double*	lowBound = problem->lowBound;
	double*	uppBound = problem->uppBound;

	if (n != numVar) {
		fprintf (stderr, "Initializing problem is wrong(numVar is uncorrect)\n");
		exit (-1);
	}


	for (i=numObj-1; i<n; i++) {
		// base on paper
		// x[i] = (1+(i+2.0-numObj)/(n-numObj+1))*(x[i]-lowBound[i]) - x[0]*(uppBound[i]-lowBound[i]);

		// base on PlatEMO
		x[i] = (1+(i+1.0)/n)*(x[i]-lowBound[i]) - x[0]*(uppBound[i]-lowBound[i]);
	}
}

void L2 (double *x, int n, Problem_t *problem) {
	int 	i=0;
	int 	numObj = problem->numObj;
	int 	numVar = problem->numVar;
	double*	lowBound = problem->lowBound;
	double*	uppBound = problem->uppBound;

	if (n != numVar) {
		fprintf (stderr, "Initializing problem is wrong(numVar is uncorrect)\n");
		exit (-1);
	}

	for (i=numObj-1; i<n; i++) {
		// base on paper
		// x[i] = (1+cos(0.5*PI*(i+2.0-numObj)/(n-numObj+1)))*(x[i]-lowBound[i]) - x[0]*(uppBound[i]-lowBound[i]);

		// base on PlatEMO
		x[i] = (1+cos(0.5*PI*(i+1.0)/n))*(x[i]-lowBound[i]) - x[0]*(uppBound[i]-lowBound[i]);
	}
}
