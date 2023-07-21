#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algorithm.h"
#include "mystring.h"
#include "dominate.h"
#include "link.h"
#include "algebra.h"
#include "myrandom.h"
#include "iepsilon.h"
#include "terminal.h"
#include "recombination.h"


Population_t* dlsmoea (Problem_t *problem) {
	int nPop = Population_getSize (problem);
	Population_t *P = Population_new (problem, nPop);
	// Population_t *Q=NULL, *U=NULL; 

	double g = 0;
	double maxGen = 2.0e+4;
	int p1, p2, i;
	double q[problem->numVar];
	double f[problem->numObj];
	

	while ( !isTerminal(P)) {
		g = 0;
		while (g < maxGen) {
			for (i=0; i<nPop; i++) {
				p1 = nPop * randu ();
				p2 = nPop * randu ();
				SBX_reproduce (P->var->elements+p1*P->var->colDim, P->var->elements+p2*P->var->colDim, q, f);
				// Update P with q using the HV	indicator
					
			}
			g++;
		}
	}
	
	return NULL;
}

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
