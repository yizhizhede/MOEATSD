#include "controller.h"
#include "individual.h"
#include "dominate.h"
#include "myrandom.h"
#include "logs.h"
#include "recombination.h"
#include "string.h"
#include "link.h"
#include "parameter.h"
#include <float.h>

// get a random number lying in [0,1).
static double  ctr_next_U ();

// get a set of subgroups for production
static List_t* least_element_variable (Population20_t *pop);
static List_t* nonleast_element_variable (Population20_t *pop);
static Link_t* get_niche_variables (Population20_t *pop);

// get a set of subgrups for prediction
static List_t* get_prediction_list (Population20_t *pop);

// generate new individuals in the way of subgroups
static void control_subgroup_generate  (Population20_t *pop);
static void control_subgroup_generate  (Population20_t *pop, List_t *list, int type); 

// genereate new individuals by desire
static void control_desire_generate  (Population20_t *pop);

/************************************************************************************************************/ 
/*					Controller							    */
/************************************************************************************************************/
int Controller (Problem_t *problem, Population20_t* pop) {
	double	  r;
	int	  i, index;
	Link_t*	  link = NULL;
	List_t*	  list = NULL;
	long long tmp;

#ifndef part_release
	// 1. Progress bar
	Problem_progressbar ();
#endif

	// 2. Print population to directory './output/'
	if (Problem_isTick ()) {
		Population20_print (pop);
	}

	// 3. Meet the terminal circumstance, then return 0
	if (1 == Problem_isEnd () && pop->spare == NULL) {
		return CONTROL_RETURN;
	}

	// 4. Compute the chaos: The variable preference grouping (II)
	Population20_chaos (pop);

	// 5. exponential distribution: The Poisson process (IV) 
	Population20_exponential_distribution (pop);

	// 6. adjust the distribution index of SBX
	Population20_adjust_sbx (pop);

	// 8. control the behavior of EA
	if (pop->probability > 0.95)  {
		// 8.1 increase the size of population
		if ((r = ctr_next_U ()) < 0.1/pop->popSize && pop->popSize < 0) { 	// 512
			return CONTROL_REPRODUCE; 		
		}

		// 8.2 partial grouping using the oldest individual.
		if ((r = ctr_next_U ()) < 0.1/pop->popSize && 0 == Problem_isEnd ()) {
			index = 0;
			tmp = pop->inds[0]->id;
			for (i=1; i<pop->popSize; i++) {
				if (pop->inds[i]->id < tmp) {
					index = i;
					tmp = pop->inds[i]->id;
				}
			}
			Population20_partial (pop, pop->inds[index]);	// 8. partial order gouping (I)
			return CONTROL_REPRODUCE_ELIMINATE;
		}

		// 8.3 principal component analysis
		if ((r = ctr_next_U ()) < 1.0/pop->popSize && 0 == Problem_isEnd ()) { // 9. principal component analysis
			// get the subgroup
			list = get_prediction_list (pop);

			// select a subgroup randomly 
			i = (int)((ctr_next_U ())*list->nLink);
			link = list->list_head;
			while (i>0) {
				i--;	
				link = link->next;
			}
				
			// PCA
			Population20_principal (pop, link);		// PCA (III)
		
			//
			List_free (&list);
			return CONTROL_REPRODUCE_ELIMINATE;
		}


		// 8.4 generate new individual by desire
		if ( randu () < 1.0/pop->popSize && 0 == Problem_isEnd ()) { 
			control_desire_generate  (pop);
			return CONTROL_REPRODUCE_ELIMINATE;
		}
	}

	// 9. reproduce and eliminate
	return CONTROL_REPRODUCE_ELIMINATE;
}

/************************************************************************************************************/ 
/*					control_return							    */
/************************************************************************************************************/
void control_return (Population20_t *pop) {
	// 1. Print pop to directory ./output/
	Population20_print (pop);

	// 2. Free population
	Population20_free (&pop);
}

/************************************************************************************************************/ 
/*					control_reproduce						    */
/************************************************************************************************************/
void control_reproduce (Population20_t *pop) {
	// if the spare is empty, generate new individuals
	while (NULL == pop->spare && 0 == Problem_isEnd ()) {
		control_subgroup_generate (pop);
	}

	//
	if (pop->spare != NULL) {
		// add one single individual into population.
		Population20_add (pop, pop->spare->pInd);	

		// delede the individual from pop->spare
		Individual_node_del (&pop->spare, pop->spare->pInd);
	}
}

/************************************************************************************************************/ 
/*					control_eliminate						    */
/************************************************************************************************************/
void control_eliminate (Population20_t *pop) {
	Individual_t *ind=NULL;
	
	// Step 1. Find the worse individual in populaiton
	ind = Population20_worst (pop);
	
	// Step 2. Reduce a individual
	Population20_del (pop, ind);	

	// Step 3. Free the memory of the detached individual.
	Individual_free (&ind);
}

static List_t* least_element_variable (Population20_t *pop) {
	List_t*	list = NULL;
	Link_t*	link = NULL;	
	int 	i;
	int	numVar = pop->inds[0]->numVar;

	// 1. The least-element variable
	for (i=0; i<numVar; i++) {
		if (1 == pop->aggressive[i]) {
			Link_add (&link, i);
		}
	}
	if (link != NULL) {
		List_add (&list, link);
		Link_free (&link);
	}

	// 2. The single preference variables
	for (i=0; i<numVar; i++){
		if (0 == pop->chaos[i]) {
			Link_add (&link, i);
		}
	}

	if (link != NULL) {
		List_add (&list, link);
		Link_free (&link);
	}
	return list;
}

static List_t* nonleast_element_variable (Population20_t *pop) {
	List_t*	list = NULL;
	Link_t*	link = NULL;	
	int 	i;
	int 	numVar = pop->inds[0]->numVar;

	// 1. Non-least-element variables
	for (i=0; i<numVar; i++){
		if (-1 == pop->aggressive[i]) {
			Link_add (&link, i);
		}
	}
	if (link != NULL) {
		List_add (&list, link);
		Link_free (&link);
	}

	// 2. Multiple preference variables
	for (i=0; i<numVar; i++){
		if (1 == pop->chaos[i]) {
			Link_add (&link, i);
		}
	}

	if (link != NULL) {
		List_add (&list, link);
		Link_free (&link);
	}

	return list;
}

static Link_t* get_niche_variables (Population20_t *pop) {
	Link_t*		link = NULL;	
	int 		i;
	int		numVar = pop->inds[0]->numVar;
	Matrix_t*	M = Matrix_new (numVar, 3);
	int*		index = NULL;
	int		threshold;

	for (i=0; i<numVar; i++) {
		if (1 == pop->aggressive[i]) {
			M->elements[i*3+0] = 1.0;
			M->elements[i*3+1] = pop->entropy[i];
		} else {
			M->elements[i*3+0] = 0;
			M->elements[i*3+1] = -pop->entropy[i];
		}
		M->elements[i*3+2] = randu ();
	}
	index = sort (M);
	threshold = (numVar < 20) ?  numVar : 20; 	// 10, 20, 100
	for (i=0; i<threshold; i++) {
		Link_add (&link, index[(numVar-1)-i]);
	}
	
	free (index);
	Matrix_free (&M);

	return link;
}

static void control_subgroup_generate  (Population20_t *pop) {
	List_t*	list = NULL;
	Link_t* link = NULL;
	int	numVar = pop->inds[0]->numVar;
	int 	numLEV = 0;
	int	i;

	for (i=0; i<numVar; i++) if (1 == pop->aggressive[i]) {
			numLEV++;
	}

	if (randu () < 1.0*numLEV / numVar) {				
		// 
		list = least_element_variable (pop);				// least-element-variables 	
		if (list != NULL) {
			if (randu () < 0.5 )
				control_subgroup_generate  (pop, list, 1);	// by time + least-element variables
			else
				control_subgroup_generate  (pop, list, 3);	// by hv + least-element variables
			List_free (&list);
		} else {
			list = nonleast_element_variable (pop);		// non-least element variables	
			if (randu () < 0.5 )
				control_subgroup_generate  (pop, list, 1);	// by time + nonleast-element variables
			else
				control_subgroup_generate  (pop, list, 3);	// by hv + nonleast-element variables
			List_free (&list);
		}
	} else {							
		// 
		list = nonleast_element_variable (pop);			// non-least-element	
		if (list != NULL) {
			if (randu () < 0.5 )
				control_subgroup_generate  (pop, list, 1);	// by time + nonleast-element variables
			else
				control_subgroup_generate  (pop, list, 3);	// by hv + nonleast-element variables
			List_free (&list);
		} else {
			list = least_element_variable (pop);		// least-element variables	
			if (randu () < 0.5 )
				control_subgroup_generate  (pop, list, 1);	// by time + least-element variables
			else
				control_subgroup_generate  (pop, list, 3);	// by hv + least-element variables
			List_free (&list);
		}
	}

	// According to niche variables, generate individuals
	link = get_niche_variables (pop);
	List_add (&list, link);						// niche variables	
	if (list != NULL) {
		if (randu () < 0.5 ) {
			control_subgroup_generate  (pop, list, 1);	// by time + niche variables
		} else {
			control_subgroup_generate  (pop, list, 3);	// by hv + niche variables
		}
		List_free (&list);
		Link_free (&link);
	} 

}

static void control_subgroup_generate  (Population20_t *pop, List_t *list, int type) {
	int 		p1, p2;
	Individual_t*	ind=NULL;
	Link_t*		link = NULL;
	Node_t*		pNode = NULL;
	int 		numVar = pop->inds[0]->numVar;
	int 		numObj = pop->inds[0]->numObj;
	int		i, j; 
	int		flag = 0;

	double 		V1[numVar+10];
	double 		V2[numVar+10];
	double 		C1[numVar+10];
	double 		C2[numVar+10];
	double 		LB[numVar+10];
	double 		UB[numVar+10];
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		len=0;
	double 		tmp=0;

	if (NULL == list) return;

	link = list->list_head;
	while (link) {
		// 1. Select two individuals from the populaiton to reproduce a new one.
		if (1 == type){
			p1 = Population20_recommand_by_time(pop);
			p2 = Population20_recommand_by_time(pop);
		} else if (2 == type) {
			p1 = Population20_recommand_by_angle(pop);
			p2 = Population20_recommand_by_angle(pop);
		} else if (3 == type) {
			p1 = Population20_recommand_by_exchv(pop);
			p2 = Population20_recommand_by_exchv(pop);
		} else {
			p1 = Population20_recommand_by_rand(pop);
			p2 = Population20_recommand_by_rand(pop);
		}

		// 2. extract a sub group of variabls 
		len = 0;
		pNode=link->link_head;
		while (pNode) {
			V1[len] = pop->inds[p1]->var[pNode->data];
			V2[len] = pop->inds[p2]->var[pNode->data];
			LB[len] = lowBound[pNode->data];
			UB[len] = uppBound[pNode->data];
			len++;

			// next node
			pNode = pNode->next;
		}
	
		// 3. SBX
		realbinarycrossover(V1, V2, C1, C2, 1.0, len, LB, UB);	
		realmutation(C1, 1.0 / len, len, LB, UB);

		// 4. check if it is changed from parents in the aspect of variable
		for (i=0; i<len; i++) {
			tmp = V1[i] - C1[i];
			if (tmp > DBL_EPSILON || tmp < -DBL_EPSILON) {
				break;
			}
		}
		if (i >= len) {
			link = link->next;	// next node
			continue;
		}

		// 5. rebuild a new individual
		ind = Individual_new ();
		memcpy (ind->var, pop->inds[p1]->var, numVar*sizeof (double));

		len = 0;
		pNode=link->link_head;
		while (pNode) {
			ind->var[pNode->data] = C1[len];
			len++;
			pNode = pNode->next;
		}

		// 6. evaluate the new individual
		Individual_evaluate (ind);

		// 7. if parents dominate the child, then give up it.
		if (isDominate (pop->inds[p1]->obj, ind->obj, ind->numObj) == 1) {
			Individual_free (&ind);
			// next link
			link = link->next;
			continue;
		}
		if (isDominate (pop->inds[p2]->obj, ind->obj, ind->numObj) == 1) {
			Individual_free (&ind);
			// next link
			link = link->next;
			continue;
		}

		// 8. check the new individual if make enouph difference in objective after evaluating
		for (i=0, flag=1; i<pop->popSize; i++) {
			for (j=0; j<numObj; j++) {
				tmp = pop->inds[i]->obj[j] - ind->obj[j];
				if (tmp > DBL_EPSILON || tmp < -DBL_EPSILON) {
					break;
				}
			}
			if (j >= numObj) {
				flag = 0;
				Individual_free (&ind);
				break;
			}
		}
		// if the new individual doesn't make enough difference in objective function 
		if (0 == flag) {
			link = link->next;	// next node
			continue;
		}

		// 9. add individual to the spare
		Individual_node_add (&pop->spare, ind);
		ind = NULL;

		// 10. next link
		link = link->next;
	}
}

#define Q_CTR_LEN 1000000
static int      Q_ctr_index = Q_CTR_LEN;
static double   Q_ctr[Q_CTR_LEN+10];
static double   ctr_next_U () {
	int     i;
       	double  U;

	//
       	if (Q_ctr_index >= Q_CTR_LEN) {
                Q_ctr_index = 0;
                srand (time (NULL));
                for (i=0; i<Q_CTR_LEN; i++) {
                         Q_ctr[i] = 1.0 * rand () / (RAND_MAX + 1.0);
                }
        }
 
        U = Q_ctr[Q_ctr_index++];
        return U;
}

static List_t* get_prediction_list (Population20_t *pop) {
	List_t*	list=NULL;
	Link_t*	link = NULL;	
	int 	i, j;
	int 	numVar = pop->inds[0]->numVar;
	int	queue[numVar+10];
	int	tail=0;

	// 1. non-least-element variable
	for (i=0, tail=0; i<numVar; i++){
		if (-1 ==  pop->aggressive[i]) {
			queue[tail++] = i;
		}
	}

	// 2. least-element variable
	for (i=0; i<numVar; i++){
		if (1 == pop->aggressive[i]) {
			Link_add (&link, i);
			for (j=0; j<tail; j++) {
				Link_add (&link, queue[j]);
			}
			List_add (&list, link);
			Link_free (&link);
		}
	}

	// 3. The multiple preference variables
	for (i=0, tail=0; i<numVar; i++){
		if (1 ==  pop->chaos[i]) {
			queue[tail++] = i;
		}
	}

	// 4. The single preference variables
	for (i=0; i<numVar; i++){
		if (0 == pop->chaos[i]) {
			Link_add (&link, i);
			for (j=0; j<tail; j++) {
				Link_add (&link, queue[j]);
			}
			List_add (&list, link);
			Link_free (&link);
		}
	}

	// 5. return 
	return list;
}

// genereate new individuals by desire
static void control_desire_generate  (Population20_t *pop) {
	Individual_t*	ind=NULL;
	int 		i, j;
	int		index;
	double 		xValue;
	Matrix_t*	M = pop->frequence;
        double*         lowBound = Problem_getLowerBound ();
        double*         uppBound = Problem_getUpperBound ();

	ind = Individual_new ();
	for (i=0; i<ind->numVar; i++) {
	 	for (j=0, xValue=-1, index=0; j<M->rowDim; j++) {
			if (M->elements[j*M->colDim+i] > xValue) {
				xValue = M->elements[j*M->colDim+i];
				index = j;
			}
		}
		ind->var[i] = (index + randu ())*(uppBound[i] - lowBound[i])/M->rowDim + lowBound[i];
	}
	Individual_evaluate (ind);
	Individual_node_add (&pop->spare, ind);
}
