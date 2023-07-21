#ifndef PARAMETER_H
#define PARAMETER_H

typedef struct Parameter_tag {
	char 	  algorithm[128];	// 1. The title of the algorithm
	char 	  problem[128];		// 2. The name of problem
	int  	  numVar;		// 3. The number of decision variables
	int  	  numObj;		// 4. The number of objective functions
	int 	  popSize;		// 5. The population size
	long long lifetime;		// 6. The terminal condition

	// reproduction parameter
	int 	id_cx;			// 7. The index of crossover
	int 	id_mu;			// 8. The index of mutation

	// the number of run
	int	run;			// 9. The number of runing
} Parameter_t;

void 		Parameter_load (int argc, char **argv);
void 		Parameter_print ();
Parameter_t*	Parameter_get ();


#endif
