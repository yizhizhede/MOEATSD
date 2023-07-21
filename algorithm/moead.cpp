#include "shape.h"
#include "population.h"
#include "algebra.h"
#include "recombination.h"
#include "terminal.h"
#include "myrandom.h"
#include <string.h>

static double g_te (double *f, double *lambda, double *z, int len);

Population_t* moead (Problem_t *problem) {
	int i, j, k, l;
	Matrix_t* 	Lambda = Population_reference (problem); 
	Population_t*	pop = Population_new (problem, Lambda->rowDim);

	int* B[Lambda->rowDim];
	int  T = 10;
	Matrix_t* dist = Matrix_new (Lambda->rowDim);
	
	double z[problem->numObj];
	double c1[problem->numVar];
	double y1[problem->numObj];

	// Initialization
	for (i=0; i<problem->numObj; i++) {
		z[i] = 1.0e+10;
	}

	for (i=0; i<Lambda->rowDim; i++) {
		dist->elements[i*dist->colDim+i] = 0;
		for (j=i+1; j<Lambda->rowDim; j++) {
			dist->elements[i*dist->colDim+j] = 
				distance_p2p (Lambda->elements+i*Lambda->colDim, 
					      Lambda->elements+j*Lambda->colDim, Lambda->colDim);
			dist->elements[j*dist->colDim+i] = 
					dist->elements[i*dist->colDim+j];  
		}
	}
	for (i=0; i<Lambda->rowDim; i++) {
		B[i] = sort (dist, i);
	}
	

	// Update
	do {
		for (i=0; i<Lambda->rowDim; i++) {
			k = T * randu ();
			l = T * randu ();	
			k = B[i][k];
			l = B[i][l];
			SBX_reproduce (pop->var->elements+k*problem->numVar, 
     				   pop->var->elements+l*problem->numVar, c1, y1);

			for (j=0; j<problem->numObj; j++) {
				z[j] = y1[j] < z[j] ? y1[j] : z[j];
			}

			for (j=0; j<T; j++) {
				k = B[i][j];
				if (g_te (y1, Lambda->elements+k*Lambda->colDim, z, Lambda->colDim) <
				    g_te (pop->obj->elements+k*problem->numObj, 
					      Lambda->elements+k*Lambda->colDim, z, Lambda->colDim)) {
					memcpy (pop->obj->elements+k*problem->numObj, y1, problem->numObj*sizeof (double));
					memcpy (pop->var->elements+k*problem->numVar, c1, problem->numVar*sizeof (double));
				}
			}
			
			
		}

	} while (!isTerminal (pop));

	return pop;
}

static double g_te (double *f, double *lambda, double *z, int len) {
	double xValue, t;
	int i;

	for (i=0, xValue = 0; i<len; i++) {
		t = lambda[i]*(f[i] - z[i]);
		if (t > xValue)
			xValue = t;
	}

	return xValue;
}
