#include "individual.h"
#include "problem.h"
#include "myrandom.h"
#include "parameter.h"
#include "recombination.h"
#include "dominate.h"
#include "hv.h"
#include "igd.h"
#include "algebra.h"
#include "logs.h"
#include "kmeans.h"
#include "singular_value_decomposition.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

/** 1. The declaration of static routines about Pareto Front (PF) */
static void init_front_related_variable (Population20_t *pop);
static void update_front_related_variable_added (Population20_t *pop, Individual_t *ind);
static void update_front_related_variable_deleted (Population20_t *pop, Individual_t *ind);

/** 2. The declaration of static routines about exclusive hypervolome (HV) */
static void update_hv_related_variable (Population20_t *pop);

/** 3. The declaration of static routines for limited set */
static void limited_set_free (Individual_t *ind);
static void limited_set_init (Population20_t *pop, Individual_t *ind);
static int  limited_set_add (Individual_t *ind, Individual_t *add);
static int  limited_set_del (Individual_t *ind, Individual_t *del);
static int  is_in_history (Population20_t *pop, Individual_t *ind);
static Matrix_t* limited_set_getF1 (Individual_t *ind);
static void      limited_set_exchv (Population20_t *pop, Individual_t *ind);

/** 4. The declaration of static routines about Angle */
static void init_angle_related_variable (Population20_t *pop);
static void update_angle_related_variable_added (Population20_t *pop, Individual_t *ind);
static void update_angle_related_variable_deleted (Population20_t *pop, Individual_t *ind);

/** 5. The declaration of static routines about Frequence */
static void update_frequence_related_variable_deleted (Population20_t *pop, Individual_t *ind);

/** 6. The declaration of static routines about endpoint */
static void init_endpoint_related_variable (Population20_t *pop);
static void update_endpoint_related_variable_added (Population20_t *pop, Individual_t *ind);

/** 7. The declaration of static routines about Ash */
static void init_ash_related_variable (Population20_t *pop);
static void update_ash_related_variable_deleted (Population20_t *pop, Individual_t *ind);
static void print_ash_related_variable_deleted (Population20_t *pop);

/** 8. The declaration of static routines about Maximum and Minimum */
static void init_extreme_related_variable (Population20_t *pop);
static void update_extreme_related_variable_added (Population20_t *pop, Individual_t *ind);
static void update_extreme_related_variable_deleted (Population20_t *pop, Individual_t *ind);

/** 9. The declaration of static routines about Distance */
static void init_distance_related_variable (Population20_t *pop);
static void update_distance_related_variable_added (Population20_t *pop, Individual_t *ind);
static void update_distance_related_variable_deleted (Population20_t *pop, Individual_t *ind);

/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

/**
 1. The routines of the individual
*/
Individual_t* Individual_new () { return Individual_new (Problem_get ()); }
Individual_t* Individual_new (Problem_t *problem) {
	Individual_t *ind = NULL;

	// 1. ind
	if (NULL == (ind = (Individual_t *)calloc (1, sizeof (Individual_t)))) {
		fprintf (stderr, "%s:%d: calloc failed\n", __FILE__, __LINE__);
		exit (0);
	}

	// 2. ind->var 
	if (NULL == (ind->var = (double *)malloc (problem->numVar*sizeof (double)))) {
		fprintf (stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
		exit (0);
	}
	
	// 3. ind->numVar
	ind->numVar = problem->numVar;
	
	// 4. ind->obj 
	if (NULL == (ind->obj = (double *)malloc (problem->numObj*sizeof (double)))) {
		fprintf (stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
		exit (0);
	}

	// 5. ind->numObj
	ind->numObj = problem->numObj;

	return ind;	
}

void Individual_init (Individual_t *ind) {
	Problem_t*	problem = Problem_get ();
	int 		i;

	// 1. Initialize the value of variables randomly.
	for (i=0; i<ind->numVar; i++) {
		ind->var[i] = problem->lowBound[i] + randu()*(problem->uppBound[i]-problem->lowBound[i]);
	}

	// 2. Evaluate the individual
	Individual_evaluate (ind);
}

void Individual_evaluate (Individual_t *ind) {
	// 1. Evalute the individual
	Problem_evaluate (ind->var, ind->numVar, ind->obj, ind->numObj);

	// 2. Set the identity
	ind->id = Problem_getFitness ();
}

void Individual_free (Individual_t **ind) {
	if (NULL == ind || NULL == (*ind))
		return;
	if (NULL != (*ind)->var) 
		free ((*ind)->var);
	if (NULL != (*ind)->obj) 
		free ((*ind)->obj);
	if (NULL != (*ind)->dominates) 
		Individual_node_free (&((*ind)->dominates));
	if (NULL != (*ind)->dominated) 
		Individual_node_free (&((*ind)->dominated));
	if (NULL != (*ind)->innerInds) 
		limited_set_free (*ind);

	free (*ind);
	*ind = NULL;
}


Individual_t* Individual_crossover (Individual_t *ind1, Individual_t *ind2) {
	Individual_t *ind = Individual_new ();	
	Individual_t *ind3 = Individual_new ();	

	// 1. Get boundary
	double *low = Problem_getLowerBound ();
	double *upp = Problem_getUpperBound ();

	// 2. Get index of mutation and crossover, and set the value of them 
	double id_mu = (Parameter_get ())->id_mu;
	double id_cx = (Parameter_get ())->id_cx;
	SBX_setPara (id_mu, id_cx);

	// 3. Crossover
	realbinarycrossover(ind1->var, ind2->var, ind->var, ind3->var, 1.0, ind->numVar, low, upp);

	// 4. Free ind3
	Individual_free (&ind3);

	// 5. Return ind.
	return ind;
}

void Individual_mutation (Individual_t *ind) {
	double *low = Problem_getLowerBound ();
	double *upp = Problem_getUpperBound ();

	// 1. Mutation
 	realmutation(ind->var, 1.0/ind->numVar, ind->numVar, low, upp);
}

void Individual_show (Individual_t *ind) {
	int i;

	printf ("id: %lld\n", ind->id);
	printf ("front: %d\n", ind->front);
	printf ("obj: ");
	for (i=0; i<ind->numObj; i++){
		printf ("%f ", ind->obj[i]);
	}
	printf ("\n");

	printf ("exchv: %f\n", ind->exchv);
}

/**
 2. The routines of individual's nodes
*/
void Individual_node_add (Individual_node_t **link, Individual_t *ind) {
	Individual_node_t *node = (Individual_node_t *)malloc (sizeof (Individual_node_t));
	node->pInd = ind;
	node->next = NULL;

	if (NULL == (*link)) {
		(*link) = node;
	} else {
		node->next = (*link);
		(*link) = node;
	}
}

void Individual_node_app (Individual_node_t **link, Individual_t *ind) {
	Individual_node_t *pNode = NULL;
	Individual_node_t *node = (Individual_node_t *)malloc (sizeof (Individual_node_t));
	node->pInd = ind;
	node->next = NULL;

	if (NULL == (*link)) {
		(*link) = node;
	} else {
		pNode = (*link);
		while (NULL != pNode->next) {
			pNode = pNode->next;
		}
		pNode->next = node;
	}
}

void Individual_node_del (Individual_node_t **link, Individual_t *ind) {
	Individual_node_t *pNode = NULL;
	Individual_node_t *qNode = NULL;

	if (NULL == link || NULL == (*link)) return;
	
	pNode = (*link);
	if ((*link)->pInd == ind) {
		(*link) = pNode->next;
		free (pNode);
	} else {
		while (pNode->next && pNode->next->pInd != ind) 
			pNode=pNode->next;
		if (pNode->next->pInd == ind) {
			qNode = pNode->next;
			pNode->next = qNode->next;
			free (qNode);	
		}
	}
}

void Individual_node_free (Individual_node_t **link) {
	Individual_node_t *pNode = NULL;

	while (*link) {
		pNode = *link;
		*link = pNode->next;
		free (pNode);
	}
}
 
 
/**
 3. The routines of population
*/
Population20_t *Population20_new (Problem_t *problem) {
	int 		i;
	Population20_t*	pop = NULL;

	// 
	srand (time (NULL));

	// 0. Allocate the memory for a population in the heap.
	if (NULL == (pop = (Population20_t *)calloc (1, sizeof (Population20_t)))) {
		fprintf (stderr, "%s:%d: calloc failed\n", __FILE__, __LINE__);
		exit (0);
	}

	// 1. Set the popSize.
	pop->popSize = (Parameter_get())->popSize;

	// 2. Allocate the memory for the pointers of the address of the individauls.
	if (NULL == (pop->inds = (Individual_t **)calloc (pop->popSize, sizeof (Individual_t*)))) {
		fprintf (stderr, "%s:%d: calloc failed\n", __FILE__, __LINE__);
		exit (0);
	}

	// 3. Initialize the individuals. 
	for (i=0; i<pop->popSize; i++) {
		pop->inds[i] = Individual_new (problem);
		Individual_init (pop->inds[i]);
	}

	// 4. Set the the extreme values: pop->maxmum, pop->minimum. 
	init_extreme_related_variable (pop);

	// 5. Set the runtime 
	pop->startTime = clock ();
	Population20_time (pop);

        // 6. Performance indicator: IGF
        pop->igd = -1.0;


/************************************* The Parameters About individuals******************************/
	// 7. Set the value of front & dominate & dominated & nd
	init_front_related_variable (pop);

	// 8. Set the exclusive hypervolome (exchv).
	update_hv_related_variable (pop);
	
	// 9. Set the angle
	init_angle_related_variable (pop);

	//    Set the distance 
	init_distance_related_variable (pop);
/****************************************************************************************************/


	// 10. Set the frequence, chaos
	pop->frequence = Matrix_new (100, problem->numVar); 		// the number of segment is 100
	pop->chaos = (int *)calloc (problem->numVar, sizeof (int));
	pop->entropy = (double *)calloc (problem->numVar, sizeof (double));

	// 11. Set the partial order groupoing 
	pop->aggressive = (int *)calloc (problem->numVar, sizeof (int));
	Population20_partial (pop, pop->inds[0]);

	// 12. The endpoints
	init_endpoint_related_variable (pop);

	// 13. ash
	init_ash_related_variable (pop);
		
	// 14. Return pop
	return pop;
}

void Population20_free (Population20_t **pop) {
	int i;

	if (pop == NULL || (*pop)==NULL)
		return;

	// 1. Free records of F1
	if (NULL != (*pop)->added_to_F1)
		Individual_node_free(&((*pop)->added_to_F1));
	if (NULL != (*pop)->deleted_from_F1) 
		Individual_node_free(&((*pop)->deleted_from_F1));

	// 2. Free ash_of_time
	if (NULL != (*pop)->ash_of_time) 
		Link_free (&((*pop)->ash_of_time));

	// 3. Free the endpoints
	if (NULL != (*pop)->endpoints) 
		Matrix_free (&((*pop)->endpoints));

	// 4. Free memory of partial order group 
	if (NULL != (*pop)->aggressive)
		free ((*pop)->aggressive);

	// 5. Free frequence, chaos.
	if (NULL != (*pop)->frequence)
		Matrix_free (&((*pop)->frequence));	
	if (NULL != (*pop)->chaos)
		free ((*pop)->chaos);

	// 6. Free extreme values.
	if (NULL != (*pop)->maxmum)
		Matrix_free (&((*pop)->maxmum));
	if (NULL != (*pop)->minimum)
		Matrix_free (&((*pop)->minimum));

	// 7. Free individuals.
	for (i=0; i<(*pop)->popSize; i++) {
		Individual_free (&((*pop)->inds[i]));	
	}

	// 8. Free memory of pointer of individuals
	if (NULL != (*pop)->inds)
		free ((*pop)->inds);
	
	// 9. free memory of pop
	if (NULL != (*pop))
		free (*pop);

	// 10. set pop to NULL
	(*pop) = NULL;
}

void Population20_time (Population20_t *pop) {
	pop->endTime = clock ();
	pop->runtime = (double)(pop->endTime - pop->startTime) / CLOCKS_PER_SEC;
}

void Population20_igd (Population20_t *pop) {
	Matrix_t *Obj = NULL;

	// 1. get Obj
	Obj = Population20_getObj(pop);

	// 2. compute igd
	pop->igd = igd (Obj, Problem_get()->title);

	// 3. free Obj
	Matrix_free (&Obj);
}


static void init_front_related_variable (Population20_t *pop) {
	int 			i, j; 
	int			front;
	int			popSize = pop->popSize;
	int			numObj = pop->inds[0]->numObj;
	Individual_node_t*	pNode = NULL;

	// 1. Compare each one with the others
	for (i=0; i<popSize-1; i++) {
		for (j=i+1; j<popSize; j++) {
			if (1 == isDominate (pop->inds[i]->obj, pop->inds[j]->obj, numObj)) {
				Individual_node_add (&pop->inds[i]->dominates, pop->inds[j]);	
				Individual_node_add (&pop->inds[j]->dominated, pop->inds[i]);	
			} else if (1 == isDominate (pop->inds[j]->obj, pop->inds[i]->obj, numObj)) {
				Individual_node_add (&pop->inds[j]->dominates, pop->inds[i]);	
				Individual_node_add (&pop->inds[i]->dominated, pop->inds[j]);	
			}
		}
		pop->inds[i]->front = 1;	
	}
	pop->inds[i]->front = 1;

	// 2. Determine the front number
	for (front=1; front<=popSize; front++) {
		for (i=0; i<popSize; i++) {
			if (front == pop->inds[i]->front) {
				pNode = pop->inds[i]->dominates; 	// pop->inds[i] < pNode->pInd
				while (pNode) {
					pNode->pInd->front = front + 1;	// key code
					pNode = pNode->next;
				}
			}
		}
	}

	// 3. Record the changes of the F1 
	for (i=0; i<popSize; i++) if (1 == pop->inds[i]->front) {
		Individual_node_add (&pop->added_to_F1, pop->inds[i]);
	}
}

static void update_front_related_variable_added (Population20_t *pop, Individual_t *ind) {
	int 			i, front;
	int			popSize = pop->popSize;
	int			numObj = ind->numObj;
	Individual_t 		*pInd = NULL;	
	Individual_node_t 	*pNode = NULL;
	Individual_node_t 	*queue = NULL;

	// 1. compare the ind with the others, and the index of the ind must equal to popSize-1
	for (i=0; i<popSize-1; i++) {	
		if (isDominate (pop->inds[i]->obj, ind->obj, numObj)==1) {
			Individual_node_add (&pop->inds[i]->dominates, ind);	
			Individual_node_add (&ind->dominated, pop->inds[i]);	
		} else if (isDominate (ind->obj, pop->inds[i]->obj, numObj)==1) {
			Individual_node_add (&ind->dominates, pop->inds[i]);	
			Individual_node_add (&pop->inds[i]->dominated, ind);	
		}
	}

	// 2. find the max front of points that dominate the 'ind'
	front = 0;
	pNode = ind->dominated;
	while (pNode) {
		if (pNode->pInd->front > front) {
			front = pNode->pInd->front; 
		}
		pNode = pNode->next;
	}

	// 3. set the front of ind to (max + 1)
	ind->front = front + 1;
	
	// C. clean the history of record of changes. 
	Individual_node_free (&pop->added_to_F1);
	Individual_node_free (&pop->deleted_from_F1);

	// A. record the change of F1. add ind to added_to_F1.
	if (1 == ind->front) {
		Individual_node_add (&pop->added_to_F1, ind);
	}

	// 4. adjust ones dominated by the ind 
	Individual_node_app (&queue, ind);
	while (queue != NULL) {
		pInd = queue->pInd;
		Individual_node_del (&queue, pInd);	// pop one from queue

		pNode = pInd->dominates;		// pInd < pNode->pInd
		while (pNode) {
			if (pNode->pInd->front <= pInd->front) {
				pNode->pInd->front = pInd->front + 1;		// key code
				Individual_node_app (&queue, pNode->pInd);	// push one into queue

				if ( 2 == pNode->pInd->front) {	// when the front change from 1 to 2.
					Individual_node_add (&pop->deleted_from_F1, pNode->pInd);
				}
			}
			pNode=pNode->next;
		}
	}
}


static void update_front_related_variable_deleted (Population20_t *pop, Individual_t *ind) {
	Individual_t* 		pInd=NULL;
	Individual_node_t* 	pNode=NULL;
	Individual_node_t* 	qNode=NULL;
	Individual_node_t* 	queue=NULL;

	// 1. break up with the better ones
	pNode = ind->dominated;	// The better individuals: pNode->pInd < ind
	while (pNode) {
		Individual_node_del (&(pNode->pInd->dominates), ind);
		pNode = pNode->next;
	}

	// 2. break up with the worse ones
	pNode = ind->dominates;	// The worse individuals: ind < pNode->pInd
	while (pNode) {
		Individual_node_del (&pNode->pInd->dominated, ind);
		pNode = pNode->next;
	}

	// C. clean the history of record of changes. 
	Individual_node_free (&pop->added_to_F1);
	Individual_node_free (&pop->deleted_from_F1);

	// A. record the change of F1. delete ind from deleted_from_F1.
	if (1 == ind->front) {
		Individual_node_add (&pop->deleted_from_F1, ind);
	}
	
	// 3. adjust the relationships among the ones that worse than the ind.
	ind->front -= 1;	// decrease the number of front of the ind by one.
	Individual_node_app (&queue, ind);	// push the ind into a queue
	while (queue != NULL) {
		pInd = queue->pInd;
		Individual_node_del (&queue, pInd); 	// pop one from the queue

		pNode = pInd->dominates;	// Worse individuals: pInd < pNode->pInd
		while (pNode) {
			if (pNode->pInd->front == pInd->front + 2) {
				qNode = pNode->pInd->dominated;	// better that pNode: qNode->pInd < pNode->pInd
				while (qNode) {
					if (qNode->pInd->front ==  pInd->front + 1) {
						break;
					}
					qNode = qNode->next;
				}
				if (NULL == qNode) {
					pNode->pInd->front -= 1;	// key code
					Individual_node_app (&queue, pNode->pInd);	// push one into the queue

					if (1 == pNode->pInd->front){	 
						// B. record the changes of F1
						Individual_node_add (&pop->added_to_F1, pNode->pInd);
					}
				}
			}
			pNode=pNode->next;
		}
	}
	ind->front += 1;	// recover the value of ind's front
}

static void update_hv_related_variable (Population20_t *pop) {
	int 			i;
	int 			popSize = pop->popSize;
	Individual_node_t*	pNode=NULL;			// outside individuals
	int			flag = 0;

	// 1. Free ones delected from F1
	pNode = pop->deleted_from_F1;
	while (pNode) {
		limited_set_free (pNode->pInd);
		pNode = pNode->next;
	}

	// 2. Initialize ones added to F1
	pNode = pop->added_to_F1;
	while (pNode) {
		limited_set_init  (pop, pNode->pInd);
		limited_set_exchv (pop, pNode->pInd);
		pNode = pNode->next;
	}
	
	// 3. Update others
	for (i=0; i<popSize; i++) if (1 == pop->inds[i]->front) {
		// 3.1 Check if it is in the history
		if (is_in_history (pop, pop->inds[i])) {
			continue;
		}

		// 3.2 Check if the deleting operation has impact on pop->inds[i].
		flag = 0;
		pNode = pop->deleted_from_F1;
		while (pNode) {
			flag += limited_set_del (pop->inds[i], pNode->pInd);
			pNode = pNode->next;
		}
		if (flag) {
			limited_set_init  (pop, pop->inds[i]);
			limited_set_exchv (pop, pop->inds[i]);
		} else {
			// 3.3) Check if the adding operation has impact on pop->inds[i].
			flag = 0;
			pNode = pop->added_to_F1;
			while (pNode) {
				flag += limited_set_add (pop->inds[i], pNode->pInd);
				pNode = pNode->next;
			}
			if (flag || pop->isChangedMin + pop->isChangedMax > 0) {
				limited_set_exchv (pop, pop->inds[i]);
			} 
		}
	}
}


static void init_angle_related_variable (Population20_t *pop) {
	int 	i, j;
	int 	popSize = pop->popSize;	
	int	numObj = pop->inds[0]->numObj;
	double	angle = 0;

	// init
	for (i=0; i<popSize; i++) {
		pop->inds[i]->angle = 1.0e+100;
	}

	// update
	for (i=0; i<popSize-1; i++) {
		for (j=i+1; j<popSize; j++) {
			angle = vector_angle (pop->inds[i]->obj, pop->inds[j]->obj, numObj);	
			if (pop->inds[i]->front >= pop->inds[j]->front && angle < pop->inds[i]->angle) {
				pop->inds[i]->angle = angle;
			}
			if (pop->inds[j]->front >= pop->inds[i]->front && angle < pop->inds[j]->angle) {
				pop->inds[j]->angle = angle;
			}
		}
	}
}

static void update_angle_related_variable_added (Population20_t *pop, Individual_t *ind) {
	int 	i;
	int 	popSize = pop->popSize;	
	int	numObj = ind->numObj; 
	double	angle = 0;

	// init
	ind->angle = 1.0e+100;

	// update
	for (i=0; i<popSize-1; i++) {
		angle = vector_angle (pop->inds[i]->obj, ind->obj, numObj);	
		if (pop->inds[i]->front >= ind->front && angle < pop->inds[i]->angle) {
			pop->inds[i]->angle = angle;
		}
		if (ind->front >= pop->inds[i]->front && angle < ind->angle) {
			ind->angle = angle;
		}
	}
}

static void update_angle_related_variable_deleted (Population20_t *pop, Individual_t *ind) {
	int 	i, j;
	int 	popSize = pop->popSize;	
	int 	numObj = ind->numObj;
	double	angle = 0, d = 0;

	for (i=0; i<popSize; i++) {
		if (ind->front <= pop->inds[i]->front) {
			angle = vector_angle (pop->inds[i]->obj, ind->obj, numObj);	
			d = angle - pop->inds[i]->angle;
			if (d < DBL_EPSILON) {
				pop->inds[i]->angle = 1.0e+100;
				for (j=0; j<popSize; j++) {
					if (j == i || pop->inds[j]->front > pop->inds[i]->front)
						continue;
					angle = vector_angle (pop->inds[i]->obj, pop->inds[j]->obj, numObj);
					if (angle < pop->inds[i]->angle) {
						pop->inds[i]->angle = angle;
					}
				}
			}
		}
	}
}

void Population20_add (Population20_t *pop, Individual_t *ind) {
	// 1. pop->popSize
	pop->popSize += 1;

	// 2. pop->inds
	pop->inds = (Individual_t **)realloc (pop->inds, (pop->popSize)*sizeof (Individual_t *));
	pop->inds[pop->popSize-1] = ind;

	// 3. maxmum and minimum
	update_extreme_related_variable_added (pop, ind);

	// 4. update the variable relevant to PF
	update_front_related_variable_added (pop, ind);

	// 5. update the variable relevant to HV.
	if (1 == ind->front) {
		update_hv_related_variable (pop);
	}

	// 6. update the variables involving Angle.
	update_angle_related_variable_added (pop, ind);

	// update the variables involving Distance.
	update_distance_related_variable_added (pop, ind);

	// update endpints 
	update_endpoint_related_variable_added (pop, ind);
}

// detach the individual from the population
void Population20_del (Population20_t *pop, Individual_t *ind) {
	int 	i; 
	int 	popSize = pop->popSize;

	/* 1. update the pop->popSize & pop->inds */
	for (i=0; i<popSize; i++) {
		if (pop->inds[i] == ind) {
			pop->inds[i] = pop->inds[popSize-1];
			pop->popSize -= 1;
			break;
		}
	}

	/** 2. update the maxmim and minimum */
	update_extreme_related_variable_deleted (pop, ind);

	/** 3. update the varialbes relevant to Pareto Front */
	update_front_related_variable_deleted (pop, ind);

	/** 4. update the varialbes relevant to HV */
	if (1 == ind->front) {
		update_hv_related_variable (pop);
	}

	/** 5. update the varialbes relevant to Angle */
	update_angle_related_variable_deleted (pop, ind);

	/** update the varialbes relevant to Distance */
	update_distance_related_variable_deleted (pop, ind);

	/** 6. update the variable relevant to Frequence */
	update_frequence_related_variable_deleted (pop, ind);

	/** 7. update the variable relevant to Ash */
	update_ash_related_variable_deleted (pop, ind);
}

static int numPrint;
void Population20_print (Population20_t *pop) {
	char 		fn[1024];
	long 		t = time (NULL);
	FILE*		fp = NULL;
	Matrix_t*	M = NULL, *F1 = NULL;
	int		i, j;
	double		sum;
	Parameter_t 	*parameter = Parameter_get ();
	Problem_t   	*problem = Problem_get ();

	// increasing the number.
	numPrint += 1;

if (false) {
	/* 1. print the value of decison varible */
	sprintf (fn, "output/%so%02dv%05d_%s_var_%02d_%ld%04d",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm, parameter->run, t, numPrint);
	M = Population20_getVar(pop); // calling of Population20_getVar()
	Matrix_print (M, fn);	
	Matrix_free (&M);
}

	/* 2. printf the value of obj */
	sprintf (fn, "output/%so%02dv%05d_%s_obj_%02d_%ld%04d",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm,	parameter->run, t, numPrint);
	M = Population20_getObj(pop);
	F1 = Matrix_front (M);
	Matrix_print (F1, fn); 	// calling of Population20_getObj()
	Matrix_free (&F1);
	Matrix_free (&M);

	/* 3. printf the value of fitness */
	sprintf (fn, "output/%so%02dv%05d_%s_fitness_%02d_%ld%04d",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm,	parameter->run,	t, numPrint);
	fp = fopen (fn, "w");
	fprintf (fp, "%lld\n", problem->fitness);
	fclose (fp);

	/* 4. printf the value of runtime */
	sprintf (fn, "output/%so%02dv%05d_%s_time_%02d_%ld%04d",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm,	parameter->run,	t, numPrint);
	fp = fopen (fn, "w");
	Population20_time (pop);	// calling of Population20_time ()
	fprintf (fp, "%f\n", pop->runtime);
	fclose (fp);

	/* 5. printf the value of igd */
	sprintf (fn, "output/%so%02dv%05d_%s_igd_%02d_%ld%04d",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm, parameter->run,	t, numPrint);
	fp = fopen (fn, "w");
	Population20_igd (pop);		// calling of Population20_igd () 
	fprintf (fp, "%.16f\n", pop->igd);
	fclose (fp);

if (false) {
	/* 7. printf desire matrix */
	sprintf (fn, "output/%so%02dv%05d_%s_desire_%02d_%ld%04d",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm, parameter->run,	t, numPrint);
	M = Matrix_dup (pop->frequence);
	for (j=0; j<M->colDim; j++) {
		for (i=0, sum=0; i<M->rowDim; i++) {
			sum += M->elements[i*M->colDim +j];
		}
		for (i=0; i<M->rowDim; i++) {
			M->elements[i*M->colDim+j] /= sum;
		}
	}
	Matrix_print (M, fn); 	
	Matrix_free (&M);
}
}

Matrix_t *Population20_getVar (Population20_t *pop) {
	int numVar = pop->inds[0]->numVar;
	int popSize = pop->popSize;
	Matrix_t *var = Matrix_new (popSize, numVar);
	int i;

	for (i=0; i<popSize; i++) {
		memcpy (var->elements+i*numVar, pop->inds[i]->var, numVar * sizeof (double));
	}
	return var;
}

Matrix_t *Population20_getObj (Population20_t *pop) {
	int 		numObj = pop->inds[0]->numObj;
	int 		popSize = pop->popSize;
	Matrix_t*	obj = Matrix_new (popSize, numObj);
	int i;

	for (i=0; i<popSize; i++) {
		memcpy (obj->elements+i*numObj, pop->inds[i]->obj, numObj * sizeof (double));
	}
	return obj;
}

void Population20_show (Population20_t *pop) {
	int popSize = pop->popSize;
	int i, front, count;

	printf ("The size of the current population is: %d\n", popSize);

	for (i=0; i<popSize; i++) {
		Individual_show (pop->inds[i]);
	}

	printf ("Runtime is %f\n", pop->runtime);
	printf ("igd is %f\n", pop->igd);

	for (front=1; front<=popSize; front++) {
		for (i=0, count=0; i<popSize; i++) {
			if (pop->inds[i]->front == front) {
				printf ("%d ", i);				
				count++;
			}
		}
		if (count>0)
			printf ("\n");
		else 
			break;
	}
	
	
	//
	printf ("frequence\n");
	Matrix_print (pop->frequence);

	// aggressive
	for (i=0; i<pop->inds[0]->numVar; i++) {
		printf ("%d^th: agrresive-> %d\n", i, pop->aggressive[i]);
	}
	
}

/**
  As for the first front: exclusive hypervolume (exchv) plays roles.
  As for the higher front: Angle plays roles.
*/
Individual_t* Population20_worst (Population20_t *pop) {
	int 	front=-1;
	double	angle=1.0e+100;
	double 	exchv = 1.0e+100;
	int 	index=0;
	int 	i, j, k;
	int 	popSize = pop->popSize;
	int	numObj = pop->inds[0]->numObj;
	double 	dis, t;

	// find the worst individual
	for (i=0; i<popSize; i++) {
		// protect the endpoints
		for (k=0; k<numObj; k++) {
			for (j=0, dis=0; j<numObj; j++)	{
				t = pop->endpoints->elements[k*numObj+j] - pop->inds[i]->obj[j];
				dis += (t > 0 ? t : -t);
			}
			if (dis < DBL_EPSILON) {
				break;
			}
		}
		if (k < numObj)	{ // If it is equivalent to one of endpoints.
			continue;
		}
		
		// compare
		if (pop->inds[i]->front > front) {
			front = pop->inds[i]->front;
			angle = pop->inds[i]->angle;
			exchv = pop->inds[i]->exchv * pop->inds[i]->distance;
			index = i;
		} else if (pop->inds[i]->front == front) {
			if (1 == front)	{
				if (pop->inds[i]->exchv * pop->inds[i]->distance < exchv) {
					angle = pop->inds[i]->angle;
					exchv = pop->inds[i]->exchv * pop->inds[i]->distance;
					index = i;
				}
			} else {
				if (pop->inds[i]->angle < angle) {
					angle = pop->inds[i]->angle;
					exchv = pop->inds[i]->exchv * pop->inds[i]->distance;
					index = i;
				}
			}
		}
		
	}

	return pop->inds[index];
}

void Population20_partial (Population20_t *pop, Individual_t *ind) {
	double 		epsilon = 1.0e-4;
	int		sign = -1;	
	int 		nPer  = 4;		
	int 		i, j, k;

	Individual_t* 	queue[nPer + 10];
	int		tail=0;
	Individual_t* 	pInd=NULL;

	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double 		xValue[ind->numObj+10];
	double		tmp=0;
	int		numVar = ind->numVar;
	int		numObj = ind->numObj;
	
	for (i=0; i<numVar; i++) {
/**********************************************************************************************************
		1. Perturbate the i^th variable nPer times
**********************************************************************************************************/
		for (j=1, tail=0; j<=nPer; j++) {
			sign *= -1;

			pInd = Individual_new (); 
			memcpy (pInd->var, ind->var, numVar*sizeof (double));
			pInd->var[i] = ind->var[i] + sign*j*epsilon;	// perturbate
	
			if (pInd->var[i] < lowBound[i]) { 
				if (pInd->var[i] + epsilon > lowBound[i]) {
					pInd->var[i] = lowBound[i];
				} else {
					Individual_free (&pInd);
					continue;
				}
			}

			if (pInd->var[i] > uppBound[i]) { 
				if (pInd->var[i] - epsilon < uppBound[i]) {
					pInd->var[i] = uppBound[i];
				} else {
					Individual_free (&pInd);
					continue;
				}
			}

			Individual_evaluate (pInd);
			queue[tail++] = pInd;
		}

/**********************************************************************************************************
		2. Check if the i^th variable possess the ability of agress.
**********************************************************************************************************/
		for (j=0; j<numObj; j++) {
			xValue[j] = ind->obj[j];
		}
		for (k=0; k<tail; k++) {
			for (j=0; j<numObj; j++) {
				if (queue[k]->obj[j] < xValue[j]) {
					xValue[j] = queue[k]->obj[j];
				}
			}
		}

		queue[tail++] = ind;
		for (k=0; k<tail; k++) {		
			for (j=0, tmp=0; j<numObj; j++) {
				tmp += (queue[k]->obj[j] -  xValue[j]);
			}
			if (tmp < numObj * DBL_EPSILON ) {
				pop->aggressive[i] = 1;
				break;
			}
		}
		if (k >= tail) {
			pop->aggressive[i] = -1;
		}
		tail--;

		// free memory
		for (k=0; k<tail; k++) {
			Individual_free (&queue[k]);
		}
	}
	
	// log partial order grouping into log/partial
	log_vector ((char *)"partial_order_grouping", pop->aggressive, numVar);
}

void init_endpoint_related_variable (Population20_t *pop) {
	int 		i, j, k;
	int 		numObj = pop->inds[0]->numObj;
	int 		popSize	= pop->popSize;
	Matrix_t*	M = NULL;
	Matrix_t*	T = NULL;
	int* 		index = NULL;
	
	if (NULL == pop->endpoints) {
		pop->endpoints = Matrix_new (numObj, numObj);
	}

	M = Population20_getObj (pop);
	T = Matrix_new (popSize, numObj);
	for (k=0; k<numObj; k++) {
		for (i=0; i<popSize; i++) {
			for (j=0; j<numObj; j++) {
				T->elements[i*numObj+j] = M->elements[i*numObj + (j+k)%numObj];
			}
		}
		index = sort (T);
		memcpy (pop->endpoints->elements+k*numObj, M->elements+index[0]*numObj, numObj*sizeof (double));
		free (index);
	}
	// log endpoints information into log/endpoints
	log_matrix ((char *)"endpoints", pop->endpoints);

	Matrix_free (&T);
	Matrix_free (&M);
}

void update_endpoint_related_variable_added (Population20_t *pop, Individual_t *ind) {
	int 	i, j;
	int 	numObj = ind->numObj;
	double 	d;
	int	flag=0;	

	pop->isChangedEndpoints = 0;
	for (i=0; i<numObj; i++) {
		flag = 1;
		for (j=0; j<numObj; j++) {
			d =  pop->endpoints->elements[i*numObj+(i+j)%numObj] - ind->obj[(i+j)%numObj];
			if ( d >  DBL_EPSILON) {
				flag = 0;
				memcpy (pop->endpoints->elements+i*numObj, ind->obj, numObj*sizeof (double));	
				pop->isChangedEndpoints = 1;
				break;
			} else if (d > -DBL_EPSILON) {
				continue;		
			} else {
				break;
			} 
		}
		if (flag==0)
			break;
	}
	
	// log endpoints information into log/endpoints
	log_matrix ((char *)"endpoints", pop->endpoints);
}

static double history_of_ratio = 2.0;
void update_frequence_related_variable_deleted (Population20_t *pop, Individual_t *ind) {
	Matrix_t*	M = pop->frequence;
	int 		numVar = ind->numVar;	
	int 		i, j;
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double		fitness =  Problem_getFitness ();

	// 1. if it has arrived the next millennium, clear the frequence data.
	if ((pop->ratio_ash - history_of_ratio) > DBL_EPSILON) {
	//	memset (M->elements, 0, M->rowDim*M->colDim*sizeof (double));
	} 

	// 2. update frequence.
	for (j=0; j<numVar; j++) {
		i = floor (M->rowDim*(ind->var[j] - lowBound[j])/(uppBound[j]-lowBound[j]));
		if (i >= M->rowDim ) {
			i = M->rowDim - 1;
		}
		M->elements[i*M->colDim+j] += (double)(fitness - ind->id);
	}

	// 3. record history
	history_of_ratio = pop->ratio_ash;
}

void Population20_chaos (Population20_t *pop) {
	Matrix_t*	M = pop->frequence;
	int 		rowDim = M->rowDim;
	int		colDim = M->colDim;

	int 		i, j;
	double		N = 0.0;
	double		p = 0.0;
	double*		U = pop->entropy;		// U(X): entropy

	// 1. Compute the entropy
	for (j=0; j<colDim; j++) {
		for (i=0, N=0; i<rowDim; i++) {
			N += M->elements[i*colDim+j];
		}
		for (i=0, U[j]=0; i<rowDim; i++) {
			p = M->elements[i*colDim+j] / N;
			if (p > DBL_EPSILON) {
				U[j] += p * log (p);
			}
		}
		U[j] = 0 - U[j] / log(2);		// log2();
	}
	
	// 2. Chaos by k-means
	if (pop->chaos != NULL)
		free (pop->chaos);
	pop->chaos = kmeans (U, colDim, 1);

	// 3. Log to files
	log_matrix ((char *)"frequnce", pop->frequence);
	log_vector ((char *)"chaos", U, colDim);	
	log_vector ((char *)"chaos", pop->chaos, colDim);	
}

#define Q_PRI_LEN  1000000
static int      Q_pri_index = Q_PRI_LEN;
static double   Q_pri[Q_PRI_LEN+10];
static double 	pri_next_U () {
	int 	i;
	double	U;

	if (Q_pri_index >= Q_PRI_LEN) {
		Q_pri_index = 0;
		srand (time (NULL));
		for (i=0; i<Q_PRI_LEN; i++) {
			Q_pri[i] = 1.0 * rand () / (RAND_MAX + 1.0);
		}
	}

	U = Q_pri[Q_pri_index++];
	return U;	
} 

void Population20_principal (Population20_t *pop, Link_t *link) {
	Individual_t*	ind = NULL;

	Matrix_t*	Var = NULL;
	Matrix_t*	subV = NULL;
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double*		O = NULL;	
	Matrix_t*	U = NULL;
	Matrix_t*	D = NULL;
	Matrix_t*	V = NULL;
	Matrix_t*	A = NULL;
	Matrix_t*	B = NULL;

	int 	i, j, k, flag, K;
	Node_t*	pNode=NULL;
	int 	popSize = pop->popSize;
	int 	numVar = pop->inds[0]->numVar;
	int 	len = 0;
	double 	d, r;
	double 	UB[numVar+10];
	double 	LB[numVar+10];
	char 	buf[128];
	long long tmp;

	// get variables
	Var = Population20_getVar (pop);
	subV = Matrix_new (popSize, Var->colDim);

	// extract the sub matrix
	len=0;
	pNode = link->link_head;
	while (pNode) {
		subV->elements[len] = Var->elements[pNode->data];
		UB[len] = uppBound[pNode->data];
		LB[len] = lowBound[pNode->data];
		len++;
		pNode = pNode->next;
	}
	subV->rowDim = popSize;		// set the number of row of subV.
	subV->colDim = len;		// set the number of column of subV.
	for (i=1; i<popSize; i++) {
		j=0;
		pNode = link->link_head;
		while (pNode) {
			subV->elements[i*len+j] = Var->elements[i*numVar+pNode->data];
			j++;
			pNode = pNode->next;
		}
	}

	// mean 
	O = Matrix_mean (subV);

	sprintf (buf, "mean:");
	log_string ((char *)"principal", buf);
	log_vector ((char *)"principal", O, len);

	// translation 
	for (j=0; j<len; j++) {
		for (i=0; i<popSize; i++) {
			subV->elements[i*len+j] -= O[j];
		}
	}

	// principal component analysis
	U = Matrix_new (popSize, len);
	D = Matrix_new (popSize, len);
	V = Matrix_new (popSize, len);
	svd (subV->elements, subV->rowDim, subV->colDim, U->elements, D->elements, V->elements); // A = UDV'
	if (subV->rowDim > subV->colDim) {
		U->rowDim = subV->rowDim;
		U->colDim = subV->colDim;
		D->rowDim = subV->colDim; 
		D->colDim = 1;
		V->rowDim = subV->colDim;
		V->colDim = subV->colDim;
	} else {
		U->rowDim = subV->rowDim;
		U->colDim = subV->rowDim;
		D->rowDim = subV->rowDim; 
		D->colDim = 1;
		V->rowDim = subV->colDim;
		V->colDim = subV->rowDim;
	}

	sprintf (buf, "U:");
	log_string ((char *)"principal", buf);
	log_matrix ((char *)"principal", U);
	sprintf (buf, "D:");
	log_string ((char *)"principal", buf);
	log_matrix ((char *)"principal", D);
	sprintf (buf, "V:");
	log_string ((char *)"principal", buf);
	log_matrix ((char *)"principal", V);

	// set U to the random values 
	A = Matrix_min (U);
	B = Matrix_max (U);
	for (j=0; j<U->colDim; j++) {
		// random 
		r = pri_next_U ();

		// 
		U->elements[j] = r * (B->elements[j]-A->elements[j]) + A->elements[j];
	}
	sprintf (buf, "rand U");
	log_string ((char *)"principal", buf);
	log_vector ((char *)"principal", U->elements, U->colDim);
		
	// remove the disturbance
	for (K=D->rowDim-1; K>0; K--) {
		if (D->elements[K] / D->elements[0] > 0.01) {
			break;
		}
	}
	
	// recove the subV at first row
	for (j=0; j<len; j++) {
		subV->elements[j] = 0;
		for (k=0; k<=K; k++) {
			subV->elements[j] += (U->elements[k])*(D->elements[k])*(V->elements[j*V->colDim+k]);
		}
	}
	sprintf (buf, "recovered subV");
	log_string ((char *)"principal", buf);
	log_vector ((char *)"principal", subV->elements, len);

	// anti-translation 
	for (j=0; j<len; j++) {
		subV->elements[j] += O[j];
		if (subV->elements[j] < LB[j]) {
			subV->elements[j] = LB[j];
		}
		if (subV->elements[j] > UB[j]) {
			subV->elements[j] = UB[j];
		}
	}
	sprintf (buf, "anti-translation subV");
	log_string ((char *)"principal", buf);
	log_vector ((char *)"principal", subV->elements, len);

	// generate the new individuals
	ind = Individual_new ();

	// find the oldest individual
	tmp = pop->inds[0]->id;
	j = 0;
	for (i=1; i<popSize; i++) {
		if (pop->inds[i]->id < tmp) {
			j = i;
			tmp = pop->inds[i]->id;
		}
	}

	// copy the oldest individual
	memcpy (ind->var, pop->inds[j]->var, numVar*sizeof (double));

	flag = 0;
	j=0;
	pNode = link->link_head;
	while (pNode) {
		d = ind->var[pNode->data] - subV->elements[j];
		if (d > DBL_EPSILON && -d > DBL_EPSILON) {
			flag = 1;
		}
		ind->var[pNode->data] = subV->elements[j];
		j++;
		pNode = pNode->next;
	}
	if (0 == flag) {
		Individual_free (&ind);
		Matrix_free (&Var);
		Matrix_free (&subV);
		free (O);
		Matrix_free (&U);
		Matrix_free (&D);
		Matrix_free (&V);
		Matrix_free (&A);
		Matrix_free (&B);
		return;
	}

	// evaluate the new individual
	Individual_evaluate (ind);

	// if make difference enough to live in terms of objective
	for (j=0; j<ind->numObj; j++) {
		for (i=0, flag=0; i<popSize; i++) {
			d = pop->inds[i]->obj[j] - ind->obj[j];
			if (d < -1.0e-6 || d > 1.0e-6) {
				flag = 1;
				break;
			}
		}
		if (1 == flag)
			break;
	}
	if (0 == flag) { 	// if it does't make difference enough.
		Individual_free (&ind);
		Matrix_free (&Var);
		Matrix_free (&subV);
		free (O);
		Matrix_free (&U);
		Matrix_free (&D);
		Matrix_free (&V);
		Matrix_free (&A);
		Matrix_free (&B);
		return;
	}
		
	// add individual to the spare
	Individual_node_add (&pop->spare, ind);

	// 
	sprintf (buf, "old:");
	log_string ((char *)"principal", buf);
	log_matrix ((char *)"principal", Var);

	// free memory
	Matrix_free (&Var);
	Matrix_free (&subV);
	free (O);
	Matrix_free (&U);
	Matrix_free (&D);
	Matrix_free (&V);
	Matrix_free (&A);
	Matrix_free (&B);
}


#define Q_REC_LEN  1000000
static int      Q_rec_index=Q_REC_LEN;
static double   Q_rec[Q_REC_LEN+10];
static double	rec_next_U () {
	int 	i;
	double	U;

	if (Q_rec_index >= Q_REC_LEN) {
		Q_rec_index = 0;
		srand (time (NULL));
		for (i=0; i<Q_REC_LEN; i++) {
			Q_rec[i] = 1.0 * rand () / (RAND_MAX + 1.0);
		}
	}

	U = Q_rec[Q_rec_index++];
	return U;
}
int Population20_recommand_by_time (Population20_t *pop) {
	int 		i; 
	double		r1, r2, r;
	int 		popSize = pop->popSize;
	double		fitness = Problem_getFitness ();
	Matrix_t*	M = NULL;
	int*		index=NULL;

	// 
	M = Matrix_new (popSize, 3);	// -front, age, id

	// update the age of individuals
	for (i=0; i<popSize; i++) {
		pop->inds[i]->age = (int)(fitness - pop->inds[i]->id);		// 
		M->elements[i*3+0] = -pop->inds[i]->front;
		M->elements[i*3+1] = (double)(fitness - pop->inds[i]->id); 	//  pop->inds[i]->age;
		M->elements[i*3+2] = i;
	}

	// sort
	index = sort (M);

	r1 = rec_next_U ();
	r2 = rec_next_U ();		
	r = r1 > r2 ? r1 : r2;
	r = M->elements[index[(int)(popSize * r)]*M->colDim+(M->colDim-1)];

	Matrix_free (&M);
	free (index);
	return (int)(r);
}

int Population20_recommand_by_rand (Population20_t *pop) {
	double	r;
	int 	popSize = pop->popSize;

	r = rec_next_U ();
	r = (int)(popSize * r);
	return r;
}

int Population20_recommand_by_angle (Population20_t *pop) {
	int 		i; 
	double		r1, r2, r;
	int 		popSize = pop->popSize;
	Matrix_t*	M = NULL;
	int*		index=NULL;

	// 
	M = Matrix_new (popSize, 3);	// -front, angle, id

	// 
	for (i=0; i<popSize; i++) {
		M->elements[i*3+0] = -pop->inds[i]->front;
		M->elements[i*3+1] = pop->inds[i]->angle;
		M->elements[i*3+2] = i;
	}

	// sort
	index = sort (M);

	r1 = rec_next_U ();
	r2 = rec_next_U ();		
	r = r1 > r2 ? r1 : r2;
	r = M->elements[index[(int)(popSize * r)]*M->colDim+(M->colDim-1)];

	Matrix_free (&M);
	free (index);
	return (int)(r);
}

int Population20_recommand_by_exchv (Population20_t *pop) {
	int 		i; 
	double		r1, r2, r;
	int 		popSize = pop->popSize;
	Matrix_t*	M = NULL;
	int*		index=NULL;

	// 
	M = Matrix_new (popSize, 4);	// -front, exchv, angle, id

	// 
	for (i=0; i<popSize; i++) {
		M->elements[i*4+0] = -pop->inds[i]->front;
		M->elements[i*4+1] = pop->inds[i]->exchv * pop->inds[i]->distance;
		M->elements[i*4+2] = pop->inds[i]->angle;
		M->elements[i*4+3] = i;
	}

	// sort
	index = sort (M);

	r1 = rec_next_U ();
	r2 = rec_next_U ();		
	r = r1 > r2 ? r1 : r2;
	r = M->elements[index[(int)(popSize * r)]*M->colDim+(M->colDim-1)];

	Matrix_free (&M);
	free (index);
	return (int)(r);
}

void Population20_exponential_distribution (Population20_t *pop) {
	double 		X;
	double		t;
	int 		i;
	Matrix_t*	M = NULL;
	long long 	fitness = Problem_getFitness (); 

	// find the oldest individual, and get X
	for (i=0, X=0; i<pop->popSize; i++) {
		t = (double)(fitness - pop->inds[i]->id);
		if (t > X) {
			X = t;
		}
	}
	pop->probability = 1.0 - exp(-(pop->lambda) * X);

	M = Matrix_new (1, 3);
	M->elements[0] = X;
	M->elements[1] = 1.0/ pop->lambda;
	M->elements[2] = pop->probability;
	log_string ((char *)"exponential_distribution", (char *)("X, Theta, P(X)"));
	log_matrix ((char *)"exponential_distribution", M);
	Matrix_free (&M);
}

#define SBX_CASE 7
static int sbx_flag = -1;
void  Population20_adjust_sbx (Population20_t *pop) {                  // 
         double  plan[SBX_CASE+10] = {0, 8, 16, 20, 32, 64, 128};
         int     season;
 
         if (1 == (sbx_flag *= -1)) {
                 season = (SBX_CASE / 2 + pop->num_of_millenium) % SBX_CASE;
         } else {
                 season = (SBX_CASE-1) - ((SBX_CASE / 2 + pop->num_of_millenium) % SBX_CASE);
         }
 
         sbx_set_mutation (plan[season]);
         sbx_set_crossover (plan[season]);
}

/**
 The routines of Ash
*/
static void init_ash_related_variable (Population20_t *pop) {
	int 	  i;
	int 	  popSize = pop->popSize;	

	// 1. init the number of millinium
	pop->num_of_millenium = 0;
	pop->lambda = 0.1 / popSize;
	print_ash_related_variable_deleted (pop);

	// 2. Step into next the episode
	for (i=0; i<popSize; i++) {
		Link_add (&(pop->ash_of_time), pop->inds[i]->id);
	}

	// 3. set ratio of ash into 1.0
	pop->ratio_ash = 1.0;
}

static void update_ash_related_variable_deleted (Population20_t *pop, Individual_t *ind) {
	int 	  i;
	int 	  popSize = pop->popSize;	
	long long fitness = Problem_getFitness ();
	Node_t*	  pNode = NULL;

	// 1. Update the pop->ash_of_time
	pNode = pop->ash_of_time->link_head;
	while (pNode) {
		if (pNode->data == ind->id) {
			// 1.1 delete the forgotten
			Link_del (pop->ash_of_time, pNode);
			free (pNode);
			pNode = NULL;

			// 2. check if it is empty
			if (NULL == pop->ash_of_time || NULL == pop->ash_of_time->link_head) {
				// 2.1 count the number of millenium
				pop->num_of_millenium++;
				pop->lambda = 1.0*pop->num_of_millenium / fitness;
				print_ash_related_variable_deleted (pop);

				// 2.2 step into next the episode
				for (i=0; i<popSize; i++) {
					Link_add (&(pop->ash_of_time), pop->inds[i]->id);
				}
			}
			break;
		} else {
			pNode = pNode->next;
		}
	}
	
	// 2. Compute the ratio_ash 
	pop->ratio_ash = 1.0 * pop->ash_of_time->nNode / popSize;
}

static void print_ash_related_variable_deleted (Population20_t *pop) {
	char 		fn[1024];
	FILE*		fp = NULL;
	Parameter_t 	*parameter = Parameter_get ();
	Problem_t   	*problem = Problem_get ();
	
	sprintf (fn, "output/%so%02dv%05d_%s_ash_%02d_0000",
		problem->title, problem->numObj, problem->numVar, parameter->algorithm,	parameter->run);
	fp = fopen (fn, "a");
	fprintf (fp, "%d,%.16f\n", pop->num_of_millenium, pop->lambda);
	fclose (fp);
}

/**
 The runtines of the extreme values
*/
static void init_extreme_related_variable (Population20_t *pop) {
	Matrix_t* 	M = NULL;

	M = Population20_getObj (pop);
	pop->maxmum = Matrix_max (M);
	pop->minimum = Matrix_min (M);

	// free M
	Matrix_free (&M);
}

static void update_extreme_related_variable_added (Population20_t *pop, Individual_t *ind) {
	int 	i;
	int	numObj = ind->numObj;

	pop->isChangedMax = 0;
	pop->isChangedMin = 0;
	for (i=0; i<numObj; i++) {
		if (ind->obj[i] > pop->maxmum->elements[i]) {
			pop->maxmum->elements[i] = ind->obj[i];
			pop->isChangedMax = 1;
		}
		if (ind->obj[i] < pop->minimum->elements[i]) {
			pop->minimum->elements[i] = ind->obj[i];
			pop->isChangedMin = 1;
		}
	}
}

static void update_extreme_related_variable_deleted (Population20_t *pop, Individual_t *ind) {
	int 	i, j;
	int 	popSize = pop->popSize;
	int 	numObj = ind->numObj;
	double	xValue;

	pop->isChangedMax = 0;
	pop->isChangedMin = 0;
	for (i=0; i<numObj; i++) {
		if (pop->maxmum->elements[i] - ind->obj[i] < DBL_EPSILON ) {
			pop->isChangedMax = 1;
			for (j=1, xValue = pop->inds[0]->obj[i]; j<popSize; j++) {
				if (pop->inds[j]->obj[i] > xValue) {
					xValue =  pop->inds[j]->obj[i];  
				}
			}
			pop->maxmum->elements[i] = xValue;
		}
		if (ind->obj[i] - pop->minimum->elements[i] < DBL_EPSILON) {
			pop->isChangedMin = 1;
			for (j=1, xValue = pop->inds[0]->obj[i]; j<popSize; j++) {
				if (pop->inds[j]->obj[i] < xValue) {
					xValue = pop->inds[j]->obj[i];  
				}
			}
			pop->minimum->elements[i] = xValue;
		}
	}
}

/**
 The runtines of the limited set
*/
static int limited_set_add (Individual_t *ind, Individual_t *add) {
	InnerInd_t*	pII = NULL;
	int 		i, j;	
	int		numObj = ind->numObj;
	int 		positive, equal;
	double		d;

	// 1.1 Check if the individual 'add' is equevalent to ind.
	if (ind->id == add->id) {
		return 0;
	}

	// 1.2 Check if 'add' has exited in the inner individuals of 'ind'
	for (i=0; i<ind->nii; i++) if (ind->innerInds[i]->origin == add->id) {
		return 0;
	}
	
	// 2.1 Construct a new inner individual. 
	pII = (InnerInd_t *) malloc (sizeof (InnerInd_t));
	pII->origin = add->id;
	for (j=0; j<numObj; j++) {
		pII->obj[j] = (add->obj[j] > ind->obj[j]) ? (add->obj[j]) : (ind->obj[j]); 
	}

	// 2.2 Check if ind->innerInds is empty
	if (NULL == ind->innerInds || ind->nii < 1) {
		ind->innerInds = (InnerInd_t **)malloc (sizeof (InnerInd_t*));
		ind->innerInds[0] = pII; 
		ind->nii = 1;
		return 1;
	}
	
	// 4. Check the domination relationships between pII and ind->innerInds[i]
	for (i=0; i<ind->nii; i++) {
		for (j=0, positive=0, equal=0; j<numObj; j++) {
			d = pII->obj[j] - ind->innerInds[i]->obj[j];
			if (d < -DBL_EPSILON) {
				positive++;	
			} else if (d < DBL_EPSILON) {
				equal++;
			}
		}
		if (0 == positive) {	// If the 'add' is dominated, quit it. 
			free (pII);
			return 0;
		} else if (positive + equal == numObj) { 	// If the 'add' dominates ind->innerInds[i], remove it.
			free (ind->innerInds[i]);
			ind->innerInds[i] = ind->innerInds[ind->nii-1];
			ind->nii--;
			i--;
		}
	}

	// 5. append 'add' into ind->innerInds.
	ind->innerInds = (InnerInd_t **)realloc (ind->innerInds, (ind->nii+1)*sizeof (InnerInd_t*));
	ind->innerInds[ind->nii] = pII;
	ind->nii++; 
	return 1;
}

static int limited_set_del (Individual_t *ind, Individual_t *del) {
	int 	i;	

	for (i=0; i<ind->nii; i++) {
		if (ind->innerInds[i]->origin == del->id) {
			free (ind->innerInds[i]);
			ind->innerInds[i] = ind->innerInds[ind->nii-1];
			ind->nii--;
			return 1;
		}
	}
	return 0;
}

static void limited_set_free (Individual_t *ind) {
	int 	i;	

	for (i=0; i<ind->nii; i++) {
		free (ind->innerInds[i]);
	}
	free (ind->innerInds);
	ind->innerInds = NULL;
	ind->nii = 0;
}

static void limited_set_init (Population20_t *pop, Individual_t *ind) {
	int 	i;
	int 	popSize = pop->popSize;

	for (i=0; i<popSize; i++) if (1 == pop->inds[i]->front && ind != pop->inds[i]) {
		limited_set_add (ind, pop->inds[i]);
	}
}

static int  is_in_history (Population20_t *pop, Individual_t *ind) {
	Individual_node_t*	qNode = NULL;

	// 1.
	qNode = pop->deleted_from_F1;	// 
	while (qNode) {
		if (ind == qNode->pInd) {
			return 1;
		}
		qNode = qNode->next;
	}

	// 2.
	qNode = pop->added_to_F1;	//
	while (qNode) {
		if (ind == qNode->pInd) {
			return 1;
		}
		qNode = qNode->next;
	}
	return 0;
}

static void limited_set_exchv (Population20_t *pop, Individual_t *ind) {
	Matrix_t*	M = NULL;
	Matrix_t*	N = NULL;
	double 		volume, d, inchv;
	int 		j;
	int		numObj = ind->numObj;

	// 1) get F1
	M = limited_set_getF1 (ind);

	// 2) normalize 
	N =  Matrix_norm (M, pop->maxmum, pop->minimum); // normalize 

	// 3) comput hv
	volume = hv (N);

	// 4) free M and N
	Matrix_free (&M);
	Matrix_free (&N);

	// 5) compute the inchv
	for (j=0, inchv = 1; j<numObj; j++) {
		d = pop->maxmum->elements[j] - pop->minimum->elements[j];
		if (d > DBL_EPSILON) {
			inchv *= 1.1 - (ind->obj[j] - pop->minimum->elements[j]) / d;
		} else {
			inchv *= 1.1;
		} 
	}

	// 6) compute the exchv 
	ind->exchv = inchv - volume;
}

static Matrix_t*  limited_set_getF1 (Individual_t *ind) {
	InnerInd_t**	ppII = NULL;	
	int		nii = ind->nii;
	int		numObj = ind->numObj;
	Matrix_t*	M = NULL;
	int		i, j;

	if (ind->nii < 1) {
		return NULL;
	}

	M = Matrix_new (nii, numObj);
	ppII = ind->innerInds;
	for (i=0; i<nii; i++) {
		for (j=0; j<numObj; j++) {
			M->elements[i*numObj+j] = ppII[i]->obj[j];
		}
	}
	return M;
}


/** 9. The declaration of static routines about Distance */
static void init_distance_related_variable (Population20_t *pop) {
	int 	i, j, k;
	int 	popSize = pop->popSize;	
	int	numObj = pop->inds[0]->numObj;
	double	d, t;
	double  s[numObj+10];

	// init
	for (i=0; i<popSize; i++) {
		pop->inds[i]->distance = 1.0e+100;
	}

	for (j=0; j<numObj; j++) {
		s[j] = pop->maxmum->elements[j] - pop->minimum->elements[j];
	}

	// update
	for (i=0; i<popSize; i++) if (1 == pop->inds[i]->front){
		for (j=i+1; j<popSize; j++) if (1==pop->inds[j]->front){
			for (k=0, d=0; k<numObj; k++) {
				t = (pop->inds[i]->obj[k] - pop->inds[j]->obj[k]) / (s[k]>DBL_EPSILON ? s[k] : 1);
				d += t*t;
			}
			d = sqrt (d);
			if (d < pop->inds[i]->distance){
				pop->inds[i]->distance = d;
			}
			if (d < pop->inds[j]->distance){
				pop->inds[j]->distance = d;
			}
		}
	}
}

static void update_distance_related_variable_added (Population20_t *pop, Individual_t *ind) {
	int 	i, j, k;
	int 	popSize = pop->popSize;	
	int	numObj = pop->inds[0]->numObj;
	double	d, t;
	double 	s[numObj+10];

	// init
	ind->distance = 1.0e+100;

	if (ind->front > 1) 
		return;

	if (pop->isChangedMax + pop->isChangedMin > 0) {
		init_distance_related_variable (pop);
		return;
	}

	for (j=0; j<numObj; j++) {
		s[j] = pop->maxmum->elements[j] - pop->minimum->elements[j];
	}

	// update 
	for (i=0; i<popSize-1; i++) if (1 == pop->inds[i]->front){
		for (k=0, d=0; k<numObj; k++) {
			t = (pop->inds[i]->obj[k] - ind->obj[k]) / (s[k]>DBL_EPSILON ? s[k]: 1);
			d += t*t;
		}
		d = sqrt (d);
		if (d < pop->inds[i]->distance){
			pop->inds[i]->distance = d;
		}
		if (d < ind->distance){
			ind->distance = d;
		}
	}
}

static void update_distance_related_variable_deleted (Population20_t *pop, Individual_t *ind) {
	int 	i, j, k;
	int 	popSize = pop->popSize;	
	int	numObj = pop->inds[0]->numObj;
	double	d, t;
	double 	s[numObj+10];

	if (ind->front > 1) 
		return;

	if (pop->isChangedMax + pop->isChangedMin > 0) {
		init_distance_related_variable (pop);
		return;
	}

	for (j=0; j<numObj; j++) {
		s[j] = pop->maxmum->elements[j] - pop->minimum->elements[j];
	}

	for (i=0; i<popSize; i++) if (1 == pop->inds[i]->front){
		for (k=0, d=0; k<numObj; k++) {
			t = (pop->inds[i]->obj[k] - ind->obj[k]) / (s[k]>DBL_EPSILON ? s[k]: 1);
			d += t*t;
		}
		d = sqrt (d);
		if (d - pop->inds[i]->distance < DBL_EPSILON){
			pop->inds[i]->distance = 1.0e+100;
			for (j=0; j<popSize; j++) if (1 == pop->inds[j]->front && i != j){
				for (k=0, d=0; k<numObj; k++) {
					t = (pop->inds[i]->obj[k] - pop->inds[j]->obj[k]) / (s[k]>DBL_EPSILON ? s[k]: 1);
					d += t*t;
				}
				d = sqrt (d);
				if (d < pop->inds[i]->distance){
					pop->inds[i]->distance = d;
				}
			}
		}
	}
}
