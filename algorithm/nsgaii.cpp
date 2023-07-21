#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algorithm.h"
#include "mystring.h"
#include "rank.h"
#include "terminal.h"
#include "reproduce.h"
#include "reduce.h"
#include "parameter.h"

Population_t* nsgaii (Problem_t *problem) {
//	int nPop = Population_getSize(problem);
	int nPop = (Parameter_get())->popSize;

	Population_t *P=NULL, *Q=NULL;
	Population_t *pop = Population_new (problem, nPop);

	// Pop_0
	P = rank_by_crowding (pop);	
	Population_free (&pop);

	while (!isTerminal (P)) {
		Q = Reproduce_equ_by_tournament(P); 
		Population_free (&P);

		P = Reduce_hal_by_cd (Q);	// reduce half of individual by crowding distance
		Population_free (&Q);
	}

	return P;
}
