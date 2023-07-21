#include "population.h"
#include "recombination.h"
#include "problem.h"
#include "snapshot.h"

//
static double progress_bar = 0;

//
int isTerminal (Population_t *pop) {
 	double nfitness = SBX_getNfitness ();
	Problem_t *problem = Problem_get ();

	if (nfitness > progress_bar*problem->lifetime && nfitness < problem->lifetime) {
		snapshot_click (pop);
		progress_bar += 0.05;
	}
	
	if (nfitness >= problem->lifetime) {
		snapshot_click (pop);
		progress_bar += 0.05;
		return 1;
	}

	return 0;
}
