#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algorithm.h"
#include "nsga3.h"
#include "ibea.h"
#include "dlsmoea.h"
#include "smsemoa.h"
#include "nsgaii.h"
#include "moead.h"
#include "two_arch2.h"
#include "thetadea.h"
#include "moeatsd.h"
#include "snapshot.h"

Population_t *moea (char *title, Problem_t *problem) {
	Population_t *pop = NULL;
	
	// snapshot 1
	snapshot_init (title, NULL);

	if (!strcmp ("NSGAII", title) || !strcmp ("NSGA-II", title))  	// NSGA-II
		pop =  nsgaii (problem);
	else if (!strcmp ("MOEAD", title)) 				// MOEA/D
		pop =  moead (problem);
	else if (!strcmp ("NSGA3", title) || !strcmp ("NSGAIII", title)) 	// MOEA-III
		pop =  nsga3 (problem);	
	else if (!strcmp ("TWOARCH2", title)) 				// two-arch2 
		pop =  two_arch2 (problem);	
	else if (!strcmp ("IBEA", title)) 				// IBEA
		pop =  ibea (problem);	
	else if (!strcmp ("THETADEA", title)) 				// theta-dea
		pop =  thetadea (problem);	
	else if (!strcmp ("DLSMOEA", title)) 				// DLS-MOEA
		pop =  dlsmoea (problem);	
	else if (!strcmp ("SMSEMOA", title)) 				// SMS-EMOA
		pop =  smsemoa (problem);	
	else if (!strcmp ("MOEATSD", title)) 				// MOEATSD
		pop =  moeatsd (problem);	
	else {
fprintf (stderr, "The algorithm %s haven't been completed now, and the ones being completed are as follows:\n", title);
		printf ("====================================\n");
		printf ("No.\t Algorithm \t Key\n");
		printf ("------------------------------------\n");
		printf ("1.\t NSGA-II \t NSGAII\n");
		printf ("2.\t MOEA/D \t MOEAD\n");
		printf ("3.\t NSGA-III \t NSGA3\n");
		printf ("4.\t TWO-ARCH2 \t TWOARCH2\n");
		printf ("5.\t IBEA \t\t IBEA\n");
		printf ("6.\t Theta-DEA \t THETADEA\n");
		printf ("7.\t DLS-MOEA \t DLSMOEA\n");
		printf ("8.\t SMS-EMOA \t SMSEMOA\n");
		printf ("9.\t MOEATSD  \t MOEATSD\n");
		printf ("====================================\n");
		
		exit (-1);
	}
	return pop;
}
