#ifndef _POPULATION_H
#define _POPULATION_H

#include "matrix.h"
#include "problem.h"
#include "link.h"

typedef struct Difference_tag {
	int 	index;
	double 	diff;
} Difference_t;

typedef struct Population_tag {
	Matrix_t 	*var;
	Matrix_t 	*obj;

	// variable difference: VD
	Difference_t	*var_diff;
	Link_t		*basicV;
	int		*I;		// record of distance variables and position variables
	int		*R;		// record of Relationship between decision variables
	double		*CBO;		// Correlation Between Objectives
	int		*CBO_posi;	// Correlation Between Objectives
	int		*CBO_nega;	// Correlation Between Objectives
	int		*CBO_zero;	// Correlation Between Objectives
	Matrix_t	*S;		// record of Shift (drift or displacement)
	double		*s;		// record of shift (drift or displacement)
	int		*C;		// record of Color of decision variables: (0, 1, ..., M-1)  and M
	double		*A;		// record of Angle of drift (or displacement)
	int		*module;	// record of module, 0: single module, 1: multi-module
	int		*M;		// record of modal, 0: single module, 1: multi-module
	int		*O;		// is least element
	
	//		MOEA/TSS
	int		*CV;		// conflict variables
	int		*CS;		// conflict space
	int		*NCOR;		// number of the correlation between variables and objectives
	int		*NINT;		// number of the interaction between variables

	//
	int		I_flag;		// flag of I 
	int		cursor;		// cursor of population	
	
	// flag of biased
	int		*biased;	// biased

	// position variables and distance variables to be optimized
	Link_t*		PV;
	Link_t*		DV;

	// DG2
	int		*dg;
	int		*linkage;	// degree of linkage 

	// a set of groups 
	List_t		*groups;	// the result of grouping

	//
	Matrix_t	*gradient;	// for resource allocation
	int		*neighbor;	// for optimize population
} Population_t;


enum Rank_t {Crowding, Hypervolume, Sum};	// type of rank
enum Select_t {Random, Tournament};		// type of select
enum Reproduce_t {SBX, Half};			// type of reproduce
typedef struct Optimize_tag {			// parameter of optimize
	Reproduce_t	reproduce;
	Link_t*		var;	
	Select_t	select;			
	Matrix_t*	repeller;		// the new value of decision variables are far from repeller
	int		isExec_difference;	//
	Rank_t		rank;
} Optimize_t;

Population_t *Population_new (Problem_t *problem, int nPop);
Population_t *Population_dup (Population_t *pop);
Population_t *Population_sub (Population_t *pop, int *subset, int n);
Population_t *Population_sub (Population_t *pop, Link_t *link);

void Population_free (Population_t **pop);
void Population_assess (Population_t *pop);
void Population_print (Population_t *pop, char *PRO_NAME, double runtime);
void Population_cat (Population_t **pop, Population_t *end);
void Population_layer (Population_t *pop, Link_t **L1, Link_t** L2, Link_t **L3);
void Population_dg (Population_t *pop);				// dg2
void Population_update_group (Population_t *pop);		// update variable group
void Population_group (Population_t *pop);			// variable group
void Population_update_neighbor (Population_t *pop);		// update neighber 
void Population_difference (Population_t *pop);			// relevance 
void Population_exec_difference (Population_t *pop, double *x);	// revise dicision variables according to difference
void Population_exec_line (Population_t *pop, double *x);	// revise dicision variables according to line 
void Population_remove_one (Population_t *pop);			// remove one individual by hv  
void Population_optimize_next (Population_t *pop, Optimize_t* opt);	// optimize next  

Population_t *Population_reproduce (Population_t *pop, char type);
Population_t *Population_mutation (Population_t *pop);
Population_t *Population_eliminate (Population_t *pop);
Population_t *Population_eliminate (Population_t *pop, int threshold);
Population_t *Population_Front1 (Population_t *pop);
Population_t *Population_compress (Population_t *pop);

int Population_getSize (Problem_t *problem);

Matrix_t* Population_reference (Problem_t *problem);

#endif
