#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "moea.h"

int main (int argc, char **argv)
{
	// srand
	srand (time (NULL));
#ifdef  part_release 
	printf ("The release version, MOEA (moea),  with -O3.\n");
#endif

#ifdef part_debug 
	Population_t*	pop = NULL;
	Parameter_t*	parameter = NULL;
	Problem_t*	problem = NULL;

	/* 1. load parameter */
	Parameter_load (argc, argv);	
	
	/* 2. print parameter to std out */
	Parameter_print ();

	/* 3. get a parameter */
	parameter = Parameter_get ();

	/* 4. initilize problem */
	problem = Problem_new ();

	/* 5. algorithm */
	pop = moea (parameter->algorithm, problem);

	/* 6. free pop */
	if (pop != NULL)
		Population_free (&pop);
#endif

	return 0;
}
