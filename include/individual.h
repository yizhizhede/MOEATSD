#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <time.h>
#include "problem.h"
#include "matrix.h"
#include "link.h"


struct Individual_tag;		// The invidual date type
struct Individual_node_tag;
struct InnerInd_tag;		// The inner individuals designed for exclusive hypervolume
struct Population20_tag;	// Population 2.0

typedef struct Individual_tag 		Individual_t;
typedef struct Individual_node_tag 	Individual_node_t;
typedef struct InnerInd_tag 		InnerInd_t;
typedef struct Population20_tag 	Population20_t; 

 
/****************************** 1. Individual_tag ************************************************/
struct Individual_tag {
	/** 1. Identification number, if id0 - id1 < 1, then the corresponding individuals are same. 
	At the same time, set its value to the current fitness number */
	long long id;	

	double* var;	// 2. decision variable
	int 	numVar;	// 3. the length of decision variable

	double*	obj;	// 4. corresponding objective function values
	int 	numObj;	// 5. the number of objective function

	/** The variables relevant to Pareto Front */
	int 			front;		// 6. The number of Pareto Front, the value is set from 1.
	Individual_node_t*	dominates;	// 7. The individuals who are dominated by this
	Individual_node_t*	dominated;	// 8. The individuals who dominate this 

	/** The variables relevant to hyper volume */
	InnerInd_t**		innerInds;	// inner individuals
	int			nii;		// the number of inner individuals
	double 			exchv;		// 12. exclusive hyper volume, set its initial value to zero.

	/* variable under the concept of time */
	int	age;		// 13. the age of the individual

	/* angle between this and others */
	double	angle;		// 14. the angle between this and others

	/* distance between this and others */
	double	distance;	// 15. the distance between this and others
};

/****************************** 2. Individual_node_tag ************************************************/
struct Individual_node_tag {
	Individual_t* 			pInd;	// 1. The pointer to individual
	struct Individual_node_tag*	next;	// 2. The next node
};

/****************************** 3. InnerInd_tag ************************************************/
struct InnerInd_tag {
	long long origin;		// The origin come from id of an outside indivisual
	double	  obj[20];
};


/****************************** 4. Population20_tag ************************************************/
struct Population20_tag {
	// 1. popSize
	int 		popSize;

	// 2. individuals
	Individual_t** 	inds;	

	// 3. extreme values
	Matrix_t*	maxmum;
	Matrix_t*	minimum;
	int 		isChangedMax;
	int 		isChangedMin;
	
	// 4. time-comsuming
	clock_t 	startTime; 
	clock_t 	endTime;	
	double  	runtime;

	// 5. performance indicator: IGD
	double 		igd;

	// 6. frequence and distribution
	Matrix_t*	frequence;	

	/* 7. value preference grouping: single-preference variables (0) and multi-preference variables (1)*/
	int*		chaos;			// type: 0: single-preference variables; 1: multi-perference variables. 

	// The entropy of variables.
	double*		entropy;		

	/* 8. partial order Grouping: least-element variables (1) and non-least-element variables (-1) */
	int*		aggressive;	// partial order grouping: the level of aggress

	// 9. protected points: annular priority comparison method.
	Matrix_t*	endpoints;
	int		isChangedEndpoints;

	// 10. spare individual: the queue of individuals
	Individual_node_t* spare;

	/* 11. The variables with the conception of time */
	Link_t*		ash_of_time;	// the memory to the time: How many the old individuals are still alive?
	int		num_of_millenium;	// the number of millinum;
	double		ratio_ash;	

	/* 12. The parameters of exponential distribution */
	double 		lambda;
	double		probability;

	// 13. The records about changes of F1
	Individual_node_t*	added_to_F1;
	Individual_node_t*	deleted_from_F1;
};


/**
	1. The routines of individual
*/
Individual_t* 	Individual_new ();
Individual_t* 	Individual_new (Problem_t *problem);
void 		Individual_free (Individual_t **ind);
void 		Individual_init (Individual_t *ind);
void 		Individual_evaluate (Individual_t *ind);
Individual_t* 	Individual_crossover (Individual_t *ind1, Individual_t *ind2);
void 		Individual_mutation (Individual_t *ind);
void 		Individual_show (Individual_t *ind);

/** 
	2. The routines of Individual node
*/
void Individual_node_add (Individual_node_t **link, Individual_t *Ind);
void Individual_node_app (Individual_node_t **link, Individual_t *Ind);
void Individual_node_del (Individual_node_t **link, Individual_t *Ind);
void Individual_node_free(Individual_node_t **link);

/**
	3. The routines of population 2.0
*/
Population20_t*	Population20_new (Problem_t *problem);
void 		Population20_free (Population20_t **pop);
void 		Population20_time (Population20_t *pop);
void 		Population20_igd (Population20_t *pop);
void 		Population20_print (Population20_t *pop);
void 		Population20_add (Population20_t *pop, Individual_t *ind);
void 		Population20_del (Population20_t *pop, Individual_t *ind);
Matrix_t*	Population20_getVar (Population20_t *pop);
Matrix_t*	Population20_getObj (Population20_t *pop);
void 		Population20_show (Population20_t *pop);
Individual_t*	Population20_worst (Population20_t *pop);
void		Population20_partial (Population20_t *pop, Individual_t *ind);	// grouping (I)
void 		Population20_chaos (Population20_t *pop);			// grouping (II)
void 		Population20_principal (Population20_t *pop, Link_t *link);	// prediction (III)
int 		Population20_recommand_by_time (Population20_t *pop);
int 		Population20_recommand_by_rand (Population20_t *pop);
int 		Population20_recommand_by_angle (Population20_t *pop);
int 		Population20_recommand_by_exchv (Population20_t *pop);
void 		Population20_exponential_distribution (Population20_t *pop);	// control (IV)
void 		Population20_adjust_sbx (Population20_t *pop);			// 


#endif
