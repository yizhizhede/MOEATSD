#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algorithm.h"
#include "terminal.h"
#include "reproduce.h"
#include "reduce.h"


Population_t* smsemoa (Problem_t *problem) {
	int nPop = Population_getSize (problem);
	Population_t *P = Population_new (problem, nPop);
	Population_t *Q=NULL; 

	while ( !isTerminal(P)) {
		Q = Reproduce_one_by_random (P); 
		Population_free (&P);		// free memory by hand

		P = Reduce_one_by_hv (Q); 	
		Population_free (&Q);		// free memory by hand
	}
	
	return P;
}
