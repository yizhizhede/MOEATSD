#include "snapshot.h"

static snapshot_t moment;

void snapshot_init (char *alg, Population_t *pop) {
	strcpy (moment.alg, alg);
	moment.pop = pop;
	moment.startTime = clock ();
}


void snapshot_click (Population_t *pop) {
	moment.pop = pop;

	moment.endTime = clock ();
	moment.runtime = (double)(moment.endTime - moment.startTime) / CLOCKS_PER_SEC;

 	Population_print (pop, moment.alg, moment.runtime);
}

