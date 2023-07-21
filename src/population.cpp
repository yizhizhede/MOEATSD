#include "population.h"
#include "matrix.h"
#include "myrandom.h"
#include "tournament.h"
#include "recombination.h"
#include "link.h"
#include "dominate.h"
#include "shape.h"
#include "parameter.h"
#include "igd.h"
#include "algebra.h"
#include "hv.h"
#include "iforest.h"
#include "anomaly.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define MAX(a,b) ((a)>(b)?(a):(b))

Population_t * Population_new (Problem_t *problem, int nPop) {
	Population_t*	pop = NULL;
	size_t 		size = 0;
	int 		i, j;
	double		d; 
	int		numVar = problem->numVar;
	int		numObj = problem->numObj;

	if (problem == NULL || nPop <= 0)  return NULL;
	
	size = sizeof (Population_t);
	pop = (Population_t *)malloc (size);
	memset (pop, 0, size);

	pop->var = Matrix_new (3*nPop, numVar);	// allocate 3 times of popsize for the increase of populaiton
	pop->obj = Matrix_new (3*nPop, numObj);
	pop->var->rowDim = nPop;
	pop->obj->rowDim = nPop;

	for (j=numVar-1; j>=0; j--) {
		d = (problem->uppBound[j] - problem->lowBound[j]);
		for (i=nPop-1; i>=0; i--) {
			pop->var->elements[i*numVar+j] = problem->lowBound[j] + randu()*d;
		}
	}

	Population_assess (pop);
	return pop;
}

Population_t * Population_dup (Population_t *pop) {
	Population_t *dup = NULL;
	size_t size = 0;
	
	size = sizeof (Population_t);
	dup = (Population_t *)malloc (size);
	memset (dup, 0, size);

	dup->var = Matrix_dup (pop->var);
	dup->obj = Matrix_dup (pop->obj);

	return dup;
}

Population_t * Population_sub (Population_t *pop, int *subset, int n) {
	Population_t *sub = NULL;
	size_t size = 0;
	
	size = sizeof (Population_t);
	sub = (Population_t *)malloc (size);
	memset (sub, 0, size);

	sub->var = Matrix_sub (pop->var, subset, n);
	sub->obj = Matrix_sub (pop->obj, subset, n);

	return sub;
}

Population_t * Population_sub (Population_t *pop, Link_t *link) {
	int *arr = Link2Array (link);	
	Population_t *subP = Population_sub (pop, arr+1, arr[0]);

	free (arr);
	return subP;
}

void Population_free (Population_t **pop) {
	Matrix_free (&(*pop)->var);
	Matrix_free (&(*pop)->obj);
	free (*pop);

	*pop = NULL;
}

void Population_assess (Population_t *pop) {
	int 	i;
	int 	numVar = pop->var->colDim;
	int 	numObj = pop->obj->colDim;

	for (i=pop->var->rowDim-1; i>=0; i--) {
		Problem_evaluate (pop->var->elements+i*numVar, numVar, pop->obj->elements+i*numObj, numObj);	
	}
}

static int counter=-1;
void Population_print (Population_t *pop, char *PRO_NAME, double runtime) {
	int 		popSize= pop->var->rowDim;
	int 		numVar = pop->var->colDim;
	int 		numObj = pop->obj->colDim;
	Matrix_t*	Obj = Matrix_dup (pop->obj);
	Matrix_t*	Var = Matrix_dup (pop->var);
	List_t*		list = ndSort (Obj);
	int*		arr  = Link2Array (list->list_head);
	char		outputPattern[1024];
	char* 		fn = NULL;
	long  		t = time (NULL);
	FILE*		fp = NULL;
	int		run = (Parameter_get ())->run;
	int 		i, j, k, a;

	// get F1
	pop->cursor = 0;
	for (a=arr[0]; a>0; a--) {
		i = arr[a];	
		memcpy (pop->var->elements+pop->cursor*numVar, Var->elements+i*numVar, numVar*sizeof (double));	
		memcpy (pop->obj->elements+pop->cursor*numObj, Obj->elements+i*numObj, numObj*sizeof (double));	
		pop->cursor++;
	}
	pop->var->rowDim = pop->cursor;
	pop->obj->rowDim = pop->cursor;
	popSize	= pop->cursor;

	// free
	free (arr);
	List_free (&list);
	Matrix_free (&Var);
	Matrix_free (&Obj);

	// counter
	counter++;

	// Pattern of output 
	sprintf (outputPattern, "output/%so%02dv%05d_%s_TYPE_%03d_%ld%04d", 
		(Problem_get())->title, numObj, numVar, PRO_NAME, run, t, counter);

	// var
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_var_");
	fp = fopen (fn, "w");
	k = (numVar >10) ? 10 : numVar;
	for (i=0; i<popSize; i++) {
		for (j=0; j<k; j++) {
			fprintf (fp, "%.15f ", pop->var->elements[i*numVar+j]);
		}
		fprintf (fp, "\n");
	}
	fclose (fp);
	free (fn);

	// obj
	fn = strrep (outputPattern,  (char *)"_TYPE_", (char *)"_obj_");
	Matrix_print (pop->obj, fn);
	free (fn);

	// igd
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_igd_");
	fp = fopen (fn, "w");
	fprintf (fp, "%.15f\n", igd (pop->obj, (Problem_get())->title));
	fclose (fp);
	free (fn);

	// time
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_time_");
	fp = fopen (fn, "w");
	fprintf (fp, "%f\n", runtime);
	fclose (fp);
	free (fn);

	// fitness
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_fitness_");
	fp = fopen (fn, "w");
	fprintf (fp, "%lld\n", Problem_getFitness ());
	fclose (fp);
	free (fn);
}

void Population_cat (Population_t **pop, Population_t *end) {
	if (*pop == NULL) {
		*pop = Population_dup (end);
	} else {
		Matrix_cat (&(*pop)->var, end->var);
		Matrix_cat (&(*pop)->obj, end->obj);
	}
}

void Population_layer (Population_t *pop, Link_t** L1, Link_t** L2, Link_t** L3) {
	List_t *list = ndSort (pop->obj, pop->obj->rowDim/2);	
	Link_t *link = NULL;

	int n=1;
	link = list->list_head;
	while (n < list->nLink) {
		Link_cat (L1, link);
		n++;
		link=link->next;
	}

	if (n==1 && link->nNode == pop->obj->rowDim / 2) {
		Link_cat (L1, link);
	} else {
		Link_cat (L2, link);
	}
	Link_cat (L3, *L1);
	Link_cat (L3, *L2);

	List_free (&list);
}

/*
static int isLeastElementVariable (Population_t *pop, int ip, int iv, int nPer) {
	int		numVar = pop->var->colDim;
	int		numObj = pop->obj->colDim;
	double 		d;
	int 		j, k;

	Matrix_t*	X = NULL;
	Matrix_t*	Y = NULL;
	double 		L[numObj+10];
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double		base; 

	//
	nPer = nPer < 5 ? nPer + 5 : nPer;

	X = Matrix_new (nPer+10, numVar);
	Y = Matrix_new (nPer+10, numObj);
	memcpy (L, pop->obj->elements+ip*numObj, numObj*sizeof (double));

	// perturb variable 
	base = lowBound[iv] + randu()*(uppBound[iv] - lowBound[iv]);
	for (k=0; k<nPer; k++) {
		memcpy (X->elements+k*numVar, pop->var->elements+ip*numVar, numVar*sizeof (double));
		d = base + k*(uppBound[iv] - lowBound[iv])/nPer;	
		if (d > uppBound[iv]) {
			X->elements[k*numVar+iv] = d - uppBound[iv] + lowBound[iv];
		} else {
			X->elements[k*numVar+iv] = d;
		}
	        Problem_evaluate (X->elements+k*numVar, numVar, Y->elements+k*numObj, numObj);
	}
	
	// check if there is the least element
	for (k=0; k<nPer; k++) {
		for (j=0; j<numObj; j++) {
			if (Y->elements[k*numObj+j] < L[j]) {
				L[j] = Y->elements[k*numObj+j];
			}
		}
	}

	for (k=0; k<nPer; k++) {		
		for (j=0, d=0; j<numObj; j++) {
			d += (Y->elements[k*numObj+j] -  L[j]);
		}
		if (d < numObj * DBL_EPSILON ) {
			if (isDominate (Y->elements+k*numObj, pop->obj->elements+ip*numObj, numObj) == 1) {
				memcpy (pop->var->elements+ip*numVar, X->elements+k*numVar, numVar*sizeof (double));	
				memcpy (pop->obj->elements+ip*numObj, Y->elements+k*numObj, numObj*sizeof (double));	
			}

			Matrix_free (&X);
			Matrix_free (&Y);
			return 1;
		}
	}

	Matrix_free (&X);
	Matrix_free (&Y);
	return 0;
}
*/


//
static double array_max(double *array,int len);
static double gamma_func(double d);

// differential grouping
void Population_dg (Population_t *pop) {
	int 		i, j;
	int 		numVar = (Problem_get ()) -> numVar;
	int 		numObj = (Problem_get ()) -> numObj;
	double	 	f_base;
	Matrix_t*	f_hat = NULL;
	Matrix_t*	F = NULL;
	Matrix_t*	Lambda = NULL;
	int*		Theta = NULL;
	double		x1[numVar+10];
	double		x2[numVar+10];
	double		m[numVar+10];
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double		buff[100], delta1, delta2;
	double		eInf, eSup, eta0=0, eta1=0, eps;
	int		run = Parameter_get()->run;
	long  		t = time (NULL);
	char		fn[256];
	char		ALG[64];
	FILE*		fp = NULL;

	f_hat = Matrix_new (numVar, 1);
	F = Matrix_new (numVar, numVar);
	Lambda = Matrix_new (numVar, numVar); 
	if (!pop->dg) {
		Theta = (int *)malloc (numVar*numVar*sizeof (int)); 
	} else {
		Theta = pop->dg; 
	}


	for (i=0; i<numVar; i++) {
		m[i]  = lowBound[i] +randu()*(uppBound[i] - lowBound[i]);
	}

	// compute f_base
	memcpy (x1, pop->var->elements, numVar*sizeof (double));
	memcpy (buff, pop->obj->elements, numObj*sizeof (double));
	f_base = sum(buff, numObj);

	// compute f_hat
	for (i=0; i<numVar; i++) {
		memcpy (x2, x1, numVar*sizeof (double));
		x2[i] = m[i];
        	Problem_evaluate (x2, numVar, buff, numObj);
		if ( isDominate (buff, pop->obj->elements, numObj) == 1) {
			memcpy (pop->var->elements, x2, numVar*sizeof (double));
			memcpy (pop->obj->elements, buff, numObj*sizeof (double));
		}
		f_hat->elements[i] = sum(buff, numObj);
	}

	// compuate F
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) {
			memcpy (x2, x1, numVar*sizeof (double));
			x2[i] = m[i];
			x2[j] = m[j];
			Problem_evaluate (x2, numVar, buff, numObj);
			if ( isDominate (buff, pop->obj->elements, numObj) == 1) {
				memcpy (pop->var->elements, x2, numVar*sizeof (double));
				memcpy (pop->obj->elements, buff, numObj*sizeof (double));
			}
			F->elements[i*numVar+j] = sum(buff, numObj);
			F->elements[j*numVar+i] = F->elements[i*numVar+j]; 
		}
	}

	// compute Lambda: function ISM ()
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) {
			delta1 = f_hat->elements[i] - f_base;
			delta2 = F->elements[i*numVar+j] - f_hat->elements[j];
			Lambda->elements[i*numVar+j] = fabs (delta1 - delta2);
			Lambda->elements[j*numVar+i] =Lambda->elements[i*numVar+j];
		}
	}

	// compute Theta: function DSM () 
	for (i=0; i<numVar*numVar; i++) {
		Theta[i] = 100;
	}
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) {
			buff[0]	= f_base;
			buff[1] = F->elements[i*numVar+j];
			buff[2] = f_hat->elements[i];
			buff[3] = f_hat->elements[j];

			eInf = gamma_func(2.0)*MAX(buff[0]+buff[1],buff[2]+buff[3]);
			eSup = gamma_func(sqrt(numVar))*array_max (buff, 4);	

			if (Lambda->elements[i*numVar+j] < eInf) {
				Theta[i*numVar+j] = 0;
				Theta[j*numVar+i] = 0;
				eta0 += 1;
			} else if (Lambda->elements[i*numVar+j] > eSup) {
				Theta[i*numVar+j] = 1;
				Theta[j*numVar+i] = 1;
				eta1 += 1;
			}
		}
	}
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) if (Theta[i*numVar+j] > 2) {
			buff[0]	= f_base;
			buff[1] = F->elements[i*numVar+j];
			buff[2] = f_hat->elements[i];
			buff[3] = f_hat->elements[j];

			eInf = gamma_func(2.0)*MAX(buff[0]+buff[1],buff[2]+buff[3]);
			eSup = gamma_func(sqrt(numVar))*array_max (buff, 4);	
			eps  = (eta0*eInf + eta1*eSup) / (eta0 + eta1);

			if (Lambda->elements[i*numVar+j] > eps) {
				Theta[i*numVar+j] = 1;
				Theta[j*numVar+i] = 1;
			} else {
				Theta[i*numVar+j] = 0;
				Theta[j*numVar+i] = 0;
			}
		}
	}

	//
	if (!pop->linkage)
		pop->linkage = (int *)calloc (numVar, sizeof(int));
	for (i=0; i<numVar; i++) {
		Theta[i*numVar+i] = 1;
		for (j=0; j<numVar; j++) if (1 == Theta[i*numVar+j]){
			pop->linkage[i]	+= 1;
		}
	}
	pop->dg = Theta;

	// print Theta 
if (numVar < 50) {
	strcpy (ALG, Parameter_get()->algorithm);
	sprintf (fn, "output/%so%02dv%05d_%s_theta_%03d_%ld", Problem_get()->title, numObj, numVar, ALG, run, t);
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		for (j=0; j<numVar; j++) {
			fprintf (fp, "%d ", Theta[i*numVar+j]);
		}
		fprintf (fp, "\n");
	}
	fclose (fp);
}

	//
	Matrix_free (&f_hat);
	Matrix_free (&F);
	Matrix_free (&Lambda); 
}

static double array_max(double *array,int len) {
        double 	t=array[len-1];
        int 	i;

        for (i=len-2; i>=0; i--) {
                if (array[i] > t) {
                        t = array[i];
		}
        }
        return t;
}

static double gamma_func(double d) {
     double muM = FLT_EPSILON /2.0;	
//   double muM = DBL_EPSILON /2.0;
     return (d * muM)/(1 - (d * muM));
}
     
//Partial Order Grouping
/*
void Population_pog (Population_t *pop) {
	int	*I = NULL;
	int 	numVar = pop->var->colDim;
	int 	nPer = 20;
	int	n_test = 5;
	int 	ip = 0, i, j;

	I = (int *)malloc (numVar*sizeof (double));
	for (i=0; i<numVar; i++) {
		I[i] =  isLeastElementVariable (pop, ip, i, nPer);
		if (0 == I[i]) {
			for (j=0; j<n_test; j++) {
				I[i] += isLeastElementVariable (pop, ip, i, (j+1)*10+nPer);
			}
			if (I[i] > 1) {
				I[i] = 1;
			}
		}
	}
	pop->pog = I;
}
*/

static void adjust_group (List_t *list) {
	Link_t	*link=NULL, *pLink = NULL, *qLink=NULL;
	
        printf ("before:\n");
        List_show (list);

	link = list->list_head->next;
	while (link) {
		if (1 == link->nNode) {
			pLink = link;
			qLink = pLink->next;
			while (pLink && qLink) {
				if (1 == qLink->nNode) {
					pLink->next = qLink->next;
					Link_add (&link, qLink->link_head->data);
					Link_free (&qLink);
					if (link->nNode >= 64) {
						break;
					}
					qLink = pLink->next;
				} else {
					pLink = qLink;
					if (qLink) {
						qLink = qLink->next;
					}
				}
			}
		}

		link = link->next;
	}

        printf ("after:\n");
        List_show (list);
}
void Population_update_group (Population_t *pop) {
	int 		numVar = pop->var->colDim;
	int 		numObj = pop->obj->colDim;
	int 		i, j, k;
	int*		I = pop->I;
	int*		Theta = pop->dg;
	Link_t*		link = NULL;
	int		vis[numVar];
	int		queue[numVar+10], head=0, tail=0;
	int 		xValue = 0;

	// arrenge position variables
	for (i=0; i<numVar; i++) {
		printf ("the linkage of [%d] is %d\n", i, pop->linkage[i]);	
		if (0 == I[i]) {
			queue[tail++] = i;
		}
	}
	for (i=0; i<tail; i++) {
		printf ("linkage: %d-->%d\n", queue[i], pop->linkage[queue[i]]);
	}

	if (0 == tail) {
		for (i=0, xValue=-1, j=0; i<numVar; i++) {
			if (pop->linkage[i] > xValue) {
				xValue = pop->linkage[i];
				j = i;
			}
		}
		pop->I[j] = 0;
		queue[tail++] = j;
	} 
	if (tail < numObj) {
		for (i=0; i<tail; i++)	{
			Link_add (&link, i);
		}
		List_add (&pop->groups, link);
		Link_free (&link);
		link = NULL;
	} else {
		for (k=0; k<numObj-1; k++) {
			for (i=0, xValue = -1; i<tail; i++) {
				if (pop->linkage[queue[i]] > xValue) {
					xValue = pop->linkage[queue[i]];
				}
			}
			for (i=0; i<tail; i++) {
				if (pop->linkage[queue[i]] == xValue) {
					Link_add (&link, queue[i]);
					queue[i] = queue[tail-1];			
					i--;
					tail--;
				}
			}
		}
		List_add (&pop->groups, link);
		Link_free (&link);
		link = NULL;

		// 0 -> 1: The left is set to distance variable
		for (i=0; i<tail; i++) {
			I[queue[i]] = 1;
		}
	}

	// search in manner of tree, arrenge distance variable
	memset (vis, 0, numVar*sizeof (int));
	while (true) {
		tail = head = 0;
		for (i=0; i<numVar; i++) if (!vis[i] && I[i] == 1) {
			vis[i] = 1;
			queue[tail++] = i;
			Link_add (&link, i);
			break;
		}
		while (tail-head > 0) {
			i = queue[head++];
			for (j=i+1; j<numVar; j++) if (!vis[j] && I[j] == 1 && Theta[i*numVar+j] == 1) {
				vis[j] = 1;
				queue[tail++] = j;
				Link_add (&link, j);
			}
		}
		List_add (&pop->groups, link);
		Link_free (&link);
		link = NULL;

		//
		for (i=0; i<numVar; i++) if (!vis[i] && I[i] == 1) {
			break;
		}
		if (i >= numVar) {
			break;
		}
	}
	// if there is no distance variable
	if (pop->groups->nLink < 2)for (i=0; i<numVar; i+=5) {
		for (j=i-5; j<i+5; j++)	{
			Link_add (&link, (j+numVar)%numVar);	
		}
		List_add (&pop->groups, link);
		Link_free (&link);
		link = NULL;
	}

	// adjust the list, split the big group, and merge the small groups
	adjust_group (pop->groups);
}

void Population_group (Population_t *pop) {
	// DG2
	Population_dg (pop);

	// update group
	Population_update_group (pop);
}

void Population_update_neighbor (Population_t *pop) {	// update neighber 
	int 		popSize = pop->var->rowDim;
	int 		numVar  = pop->var->colDim;
	int 		i, j, k, m;
	int*		index = NULL;
	double		d;
	Matrix_t*	D = Matrix_new (popSize, popSize);
	int		queue[numVar], n=0;

	//
	if (pop->neighbor) free (pop->neighbor);
	pop->neighbor = (int *)malloc (popSize*popSize*sizeof (int));

	// find position variables
	for (i=0; i<numVar; i++) if (!pop->I[i]) {
		queue[n++] = i;
	}

	// compute distance between individuals
	for (i=0; i<popSize-1; i++) {
		for (j=i+1; j<popSize; j++) {
			for (k=0; k<n; k++)	{
				m = queue[k];
				d = pop->var->elements[i*numVar+m] - pop->var->elements[j*numVar+m];
				D->elements[i*popSize+j] += d*d;
			}
			D->elements[j*popSize+i] = D->elements[i*popSize+j];
		}
	}

	// compute neighbor
	for (i=0; i<popSize; i++) {
		index = sort (D, i);
		for (j=0; j<popSize; j++) {
			pop->neighbor[i*popSize+j] = index[j];
		}
		free (index);
	}

	//
	Matrix_free (&D);
}

static void get_score_from_A (Population_t *pop, int* score) {
	int		numVar = pop->var->colDim;
	int 		i, j, k, a, b;
	Matrix_t*	M = NULL;
	int*		index = NULL;
	double		Sim[numVar+10];

	//
	M = Matrix_new (numVar, 2);
	for (i=0; i<numVar; i++) {
		M->elements[i*2+0] = pop->A[i];
		M->elements[i*2+1] = i + 0.1;
	}
	index = sort (M, (char *)"DES");

	for (i=0; i<numVar-1; i++) {
		a = index[i];	
		b = index[i+1];
		Sim[i] = (pop->A[a]) > DBL_EPSILON ? (pop->A[b] / pop->A[a]) : 1;
	}

	j = numVar > 10 ? 10 : numVar;
	while (j <= numVar) {
		k = anomaly (Sim, j);
		if (k >=0) {
			break;
		}
		j+=3;
	}
	if (k >=0) {
		for (i=0; i<=k; i++) {
			a = index[i];
			score[a] = 0;
		}
		for (i=k+1; i<numVar; i++) {
			b = index[i];
			score[b] = 1;
		}
	} else {
		for (i=0; i<numVar; i++) {
			score[i] = 1;
		}
	}

	// free
	Matrix_free (&M); free (index);
}

static int connected_graph_vis[1000000];
static Link_t* find_connected_graph (Population_t *pop, int *I, int degree) {
	Link_t*	link = NULL;
	Node_t*	pNode = NULL;
	int	numVar = pop->var->colDim;
	int 	i, j, k;

	for (i=0; i<numVar-1; i++) if (0 == I[i]) {
		for (j=i+1; j<numVar; j++) if (0 == I[j]) {
			if (pop->R[i*numVar+j] == degree && 0 == (connected_graph_vis[i] && connected_graph_vis[j])) {
				if (NULL == link) {
					if (0 == connected_graph_vis[i]) {
						Link_add (&link, i);
						connected_graph_vis[i] = 1;
					}
					if (0 == connected_graph_vis[j]) {
						Link_add (&link, j);
						connected_graph_vis[j] = 1;
					}
					continue;
				} 
				if (0 == connected_graph_vis[i]) {
					pNode = link->link_head;
					while (pNode) {
						k = pNode->data;
						if (pop->R[k*numVar+i] != degree) {
							break;
						}
						pNode=pNode->next;
					}
					if (NULL == pNode) {
						Link_add (&link, i);
						connected_graph_vis[i] = 1;
					}
				}
				if (0 == connected_graph_vis[j]) {
					pNode = link->link_head;
					while (pNode) {
						k = pNode->data;
						if (pop->R[k*numVar+j] != degree) {
							break;
						}
						pNode=pNode->next;
					}
					if (NULL == pNode) {
						Link_add (&link, j);
						connected_graph_vis[j] = 1;
					}
				}
			}
		}
	}
	return link;
}

static void compute_final_control_attribution (Population_t *pop) {
	int	numVar = pop->var->colDim;
	int	numObj = pop->obj->colDim;
	int	scoreA[numVar+10];
	int	I[numVar+10];
	int	i, n=0, degree;
	Link_t*	link = NULL;
	int*	arr = NULL;

	// 1. get score form A
	get_score_from_A (pop, scoreA);

	// 2. get position variables
	for (i=0; i<numVar; i++) {
		if ((0 == pop->I[i] && 0 == scoreA[i]) && pop->C[i] == numObj) {
			I[i] = 0;
			n++;
		} else {
			I[i] = 1;
		}
	}
	if (n < 1) for (i=0; i<numVar; i++) {
		if ((0 == pop->I[i] || 0 == scoreA[i]) && pop->C[i] == numObj) {
			I[i] = 0;
			n++;
		} else {
			I[i] = 1;
		}
	}
	if (n < 1) for (i=0; i<numVar; i++) {
		if (pop->C[i] == numObj) {
			I[i] = 0;
			n++;
		} else {
			I[i] = 1;
		}
	}
	if (n < 1) {
		fprintf (stderr, "ERROR: %s: %d: there is not conflit on objective functions, exit...", __FILE__, __LINE__);
		exit (0);
	}


	// 3.1
	for (i=0; i<numVar; i++) {
		pop->I[i] = 1;
	}

	// 3.2
	for (degree=7, n=0; degree >0; degree--) {
		link = 	find_connected_graph (pop, I, degree);
		while (link != NULL) {
			arr = Link2Array (link);
			for (i=arr[0]; i>0; i--) {
				pop->I[arr[i]] = 0;
				n++;
			}
			free (arr);
			Link_free (&link);
			link = find_connected_graph (pop, I, degree);
		}
	}

	// 3.3
	if (n < 1) for (i=0; i<numVar; i++) {
		pop->I[i] = I[i];
	}

	// for test
	for (i=0; i<4; i++) {
		pop->I[i] = 0;
	}
	for (i=i; i<numVar; i++) {
		pop->I[i] = 1;
	}
}

static void perturb_decision_variable (Population_t *pop) {
	int		numVar = pop->var->colDim;
	int		numObj = pop->obj->colDim;
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int 		nPer=10;		
	Matrix_t	*X  = Matrix_new (1, numVar); 
	Matrix_t	*Y  = Matrix_new (nPer*numVar, numObj); 
	Matrix_t	*L1 = Matrix_new (nPer, numObj); 
	Matrix_t	*L2 = Matrix_new (nPer, numObj); 
	Matrix_t	*L3 = Matrix_new (nPer, numObj); 
	double		x[nPer+10], y[nPer+10];

	double		d = 0, b=0; 
	int		i, j, k, r, m, n;
	char		outputPattern[1024];
	int		run = Parameter_get()->run;
	long  		t  = time (NULL);
	char*		fn = NULL;
	FILE*		fp = NULL;
	int		buff[numObj+10];
	Matrix_t*	Max = NULL;
	Matrix_t*	Min = NULL;
	Matrix_t*	O = NULL;

	// 0. Pattern of output 
	sprintf (outputPattern, "output/%so%02dv%05d_%s_TYPE_%03d_%ld%04d", 
		(Problem_get())->title, numObj, numVar, (Parameter_get())->algorithm, run, t, 1);

	// 1. alloc memeory for I, R, S, s, C, A.
	if (!pop->I) {
		pop->I = (int *)calloc(numVar, sizeof (int));
	}
	if (!pop->R) {
		pop->R = (int *)calloc(numVar*numVar, sizeof (int));
	}
	if (!pop->S) {
		pop->S = Matrix_new (numVar, numObj); 
	}
	if (!pop->s) {
		pop->s = (double *)calloc(numVar, sizeof (double));
	}
	if (!pop->C) {
		pop->C = (int *)calloc(numVar, sizeof (int));
	}
	if (!pop->A) {
		pop->A = (double *)calloc(numVar, sizeof (double));
	}

	// 2. get a point randomly
	b = randu ()/(nPer+1);
	for (i=0; i<numVar;  i++) {
		pop->var->elements[i] = lowBound[i] + b*(uppBound[i] - lowBound[i]);
	}
	Problem_evaluate (pop->var->elements, numVar, pop->obj->elements, numObj);

	// 3. perturb
	for (i=0; i<numVar; i++) {
		d = (uppBound[i] - lowBound[i])/(nPer+1);
		for (k=0; k<nPer; k++) {
			r = i*nPer+k;
			memcpy (X->elements, pop->var->elements, numVar*sizeof(double));
			X->elements[i] = (1+k)*d + pop->var->elements[i]; 
			Problem_evaluate (X->elements, numVar, Y->elements+r*numObj, numObj);
		}
	}

	// 4. compute the value of I
	for (i=0; i<numVar; i++) {
		memcpy (L1->elements, Y->elements+i*nPer*numObj, nPer*numObj*sizeof(double));

		for (j=0; j<numObj; j++) {
			y[j] = pop->obj->elements[j];
		}
		for (k=0; k<nPer; k++) {
			for (j=0; j<numObj; j++) {
				if (L1->elements[k*numObj+j] < y[j]) {
					y[j] = L1->elements[k*numObj+j];
				}
			}
		}
		for (k=0; k<nPer; k++) {
			d = distance_p2p (y, L1->elements+k*numObj, numObj);
			if (d < DBL_EPSILON ) {
				pop->I[i] = 1;
				break;
			}
		}
		d = distance_p2p (y, pop->obj->elements, numObj);
		if (d < DBL_EPSILON ) {
			pop->I[i] = 1;
		}
	}

	
	// 5. compute R 
	for (i=0; i<numVar; i++) { 	
		pop->R[i*numVar+i] = 7;
		for (j=i+1; j<numVar; j++) {
			memcpy (L1->elements, Y->elements+i*nPer*numObj, nPer*numObj*sizeof(double));
			memcpy (L2->elements, Y->elements+j*nPer*numObj, nPer*numObj*sizeof(double));
			for (k=nPer*numObj-1; k>=0; k--) {
				L3->elements[k] = L1->elements[k] - L2->elements[k];		// L3 = L1 - L2
			}

			// 5.1 check if L1 = L2
			if (norm (L3->elements, nPer*numObj) < FLT_EPSILON) {
				pop->R[i*numVar+j] = 7;
				pop->R[j*numVar+i] = 7;
				continue;
			}
				
			// 5.2 check if L2 = a*L1 + b
			for (m=0, n=0; m<numObj; m++) {
				for (k=0; k<nPer; k++) {
					x[k] = L1->elements[k*numObj+m];
					y[k] = L2->elements[k*numObj+m];
				}
				n += isLinear (x, y, nPer);
			}
			if (n >= numObj) {
				pop->R[i*numVar+j] += 2;
				pop->R[j*numVar+i] += 2;
			}

			// 5.3 check if L3 = L1 - L2 = t*A + B
			for (m=0, n=0; m<numObj; m++) {
				for (k=0; k<nPer; k++) {
					x[k] = k+1;
					y[k] = L3->elements[k*numObj+m];
				}
				n += isLinear (x, y, nPer);
			}
			if (n >= numObj) {
				pop->R[i*numVar+j] += 1;
				pop->R[j*numVar+i] += 1;
			}
		}
	}

	// 6. compute S 
	for (i=0; i<numVar; i++) {
		memcpy (L1->elements, Y->elements+i*nPer*numObj, nPer*numObj*sizeof(double));
		for (j=0; j<numObj; j++) {
			pop->S->elements[i*numObj+j] += fabs (L1->elements[j] - pop->obj->elements[j]); 
			for (k=1; k<nPer; k++) {
				pop->S->elements[i*numObj+j] += fabs(L1->elements[k*numObj+j]-L1->elements[(k-1)*numObj+j]); 
			}
			pop->S->elements[i*numObj+j] /= nPer; 
		}
	}

	// 7. compute s: note that it is not commented to use 's' to group, because 's' is a fool attribution
	Max = Matrix_max (pop->S);
	Min = Matrix_new (1, numObj);
	O   = Matrix_norm (pop->S, Max, Min);
	for (i=0; i<numVar; i++) {
		pop->s[i] = norm (O->elements+i*numObj, numObj);
	}
	Matrix_free (&O);	
	Matrix_free (&Min);	
	Matrix_free (&Max);	

	// 8. compute C
	for (i=0; i<numVar; i++) {
		for (j=0; j<numObj; j++) {
			if (pop->S->elements[i*numObj+j] < DBL_EPSILON) {
				buff[j] = 0;
			} else {
				buff[j] = 1;
			}
		}
		for (j=0, n=0; j<numObj; j++) {
			n += buff[j];	
		}
		if (n==0) {
			pop->C[i] = numObj+1;
		} else if (n==1) {
			for (j=0; j<numObj; j++) if (buff[j]) {
				pop->C[i] = j;
				break;
			}
		} else {
			pop->C[i] = numObj;
		}
	}

	// 9. compute A 
	for (i=0; i<numVar; i++) {
		memcpy (L1->elements, Y->elements+i*nPer*numObj, nPer*numObj*sizeof(double));
		O = Matrix_norm (L1);
		for (k=0, pop->A[i]=0; k<nPer-1; k++) {
			pop->A[i] += vector_angle(O->elements+k*numObj, O->elements+(k+1)*numObj, numObj); 
		}
		Matrix_free (&O);	
		pop->A[i] /= (nPer-1); 
	}

	// 10. compute final control attribution
	compute_final_control_attribution (pop);

	// Y
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_Y_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		fprintf (fp, "var: %d\n", i);
		for (k=0; k<nPer; k++) {
			r = i*nPer+k;
			for (j=0; j<numObj; j++) {
				fprintf (fp, "%f ", Y->elements[r*numObj+j]);
			}
			fprintf (fp, "\n");
		}
	}
	fclose (fp);
	free (fn);

	// I 
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_I_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		fprintf (fp, "%d\n", pop->I[i]);
	}
	fclose (fp);
	free (fn);

	// R
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_R_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		for (j=0; j<numVar; j++) {
			fprintf (fp, "(%d,%d) = %d\n", i, j, pop->R[i*numVar+j]);
		}
	}
	fclose (fp);
	free (fn);

	// S
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_S_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		for (j=0; j<numObj; j++) {
			fprintf (fp, "%.16f ", pop->S->elements[i*numObj+j]);
		}
		fprintf (fp, "\n");
	}
	fclose (fp);
	free (fn);

	// s
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_s_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		fprintf (fp, "%.16f\n", pop->s[i]);
	}
	fclose (fp);
	free (fn);

	// C
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_C_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		fprintf (fp, "%d\n", pop->C[i]);
	}
	fclose (fp);
	free (fn);

	//  A
	fn = strrep (outputPattern, (char *)"_TYPE_", (char *)"_A_");
	fp = fopen (fn, "w");
	for (i=0; i<numVar; i++) {
		fprintf (fp, "%0.16f\n", pop->A[i]);
	}
	fclose (fp);
	free (fn);

	Matrix_free (&L3);
	Matrix_free (&L2);
	Matrix_free (&L1);
	Matrix_free (&X);
	Matrix_free (&Y);
}

void Population_difference (Population_t *pop) {
	int		numVar = pop->var->colDim;
	double		d = 0; 
	int		i, j;

	// 1. perburbing variables and compute I, R, S, s, C, A. 
	perturb_decision_variable (pop);

	// 2. alloc memerry 
	if (!pop->var_diff) {
		pop->var_diff = (Difference_t*)malloc(numVar*sizeof(Difference_t));
	}

	// 3. init 
	for (i=0; i<numVar; i++) {
		pop->var_diff[i].index = i;
		pop->var_diff[i].diff  = 0;
	}

	// 4. compute var_diff
	for (i=0; i<numVar; i++) if (pop->I[i] && pop->var_diff[i].index  == i) { 		// x_i -> x_j : x_j = x_i + c
		for (j=i+1; j<numVar; j++) if (pop->I[j] && pop->var_diff[j].index == j){
			// check if L1 = L2
			if (pop->R[i*numVar+j] == 7) {
				pop->var_diff[j].index = i;
				pop->var_diff[j].diff  = 0;
				continue;
			}
			// check if L1  = A*L2 + B or L3 = L1 - L2 = t*A + B
			if (pop->R[i*numVar+j] >= 1) {
				pop->var_diff[j].index = i;
				pop->var_diff[j].diff  = 10.0;
			}
		}
	}

	// 5. collect basic variables
	for (i=0; i<numVar; i++) {
		j = pop->var_diff[i].index;
		d = pop->var_diff[i].diff;
		if (j == i) {
			Link_add (&pop->basicV, i);
		} else if (d > 2.0){
			Link_add (&pop->basicV, i);
		}
	}
	printf ("basic variables: \t");
	Link_show (pop->basicV);

	// 6. set PV, DV to be optimized
	for (i=0; i<numVar; i++) if (pop->I[i] == 0) {
		Link_add (&pop->PV, i);
	} else {
		j = pop->var_diff[i].index;
		d = pop->var_diff[i].diff;
		if (j == i) {
			Link_add (&pop->DV, i);
		} else if (d > 2.0){
			Link_add (&pop->DV, i);
		}
	}

	printf ("PV: \t");
	Link_show (pop->PV);
	printf ("DV: \t");
	Link_show (pop->DV);

	// 7. check if there is position and distance variables.
	if (NULL == pop->PV) {
		fprintf (stderr, "ERROR: %s: %d: there is not position variables.\n", __FILE__, __LINE__);
		exit (0);
	}
	if (NULL == pop->DV) {
		fprintf (stderr, "ERROR: %s: %d: there is not distance variables.\n", __FILE__, __LINE__);
		exit (0);
	}
}

// revise dicision variables according to difference
void Population_exec_difference (Population_t *pop, double *x) {
	int 	numVar = pop->var->colDim;
	double*	lowBound = Problem_getLowerBound ();
	double*	uppBound = Problem_getUpperBound ();
	int	base, i;
	double 	diff, d1, d2;

	for (i=0; i<numVar; i++) if (pop->I[i]) {
		base = pop->var_diff[i].index;
		diff = pop->var_diff[i].diff;
		d1 = uppBound[i] - lowBound[i];
		d2 = uppBound[base] - lowBound[base];
		if (base != i && diff < 2.0) {
			x[i] = d1*((x[base]-lowBound[base])/d2+diff) + lowBound[i];
		}
	}
}

void Population_exec_line (Population_t *pop, double *x) {
	int 		i, j, k, m;
	int 		popSize = pop->var->rowDim;
	int 		numVar  = pop->var->colDim;
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		r1, r2;
	int*		index;
	double		d1, d2, d3, d4;
	int		Q1[numVar], n1=0;	// position varialbes
	int		Q2[numVar], n2=0;	// distance variables

	for (i=0; i<numVar; i++) if (!pop->I[i]){
		Q1[n1++] = i;	
	} else {
		Q2[n2++] = i;
		x[i] = 0;
	}

	// step 2. 
	for (k=0; k<n1; k++) {
		j = Q1[k];
		index = sort (pop->var, j);	
		for (i=0; i<popSize; i++) if (pop->var->elements[index[i]*numVar+j] > x[j]) {
			break;
		}
		if (0 == i) {
			r1 = index[0];
			r2 = index[1];
		} else if (i >= popSize) {
			r1 = index[i-2];
			r2 = index[i-1];
		} else {
			r1 = index[i-1];
			r2 = index[i];
		}
		for (i=0; i<n2; i++) {
			m = Q2[i];
			d1 = pop->var->elements[r2*numVar+j] - pop->var->elements[r1*numVar+j];
			d2 = pop->var->elements[r2*numVar+m] - pop->var->elements[r1*numVar+m];
			d3 = x[j] - pop->var->elements[r1*numVar+j];
			d4 = pop->var->elements[r2*numVar+m] + pop->var->elements[r1*numVar+m];
			if (d1 > DBL_EPSILON || d1 < -DBL_EPSILON) {
				x[m] += d3*d2/d1+pop->var->elements[r1*numVar+m]; 
			} else {
				x[m] += d4/2; 
			}
		}
		free (index);
	}

	// step 3.
	for (i=0; i<n2; i++) {
		m = Q2[i];
		x[m] /= n1; 
		if (x[m] > uppBound[m]) {
			x[m] = uppBound[m];
		} else if (x[m] < lowBound[m]) {
			x[m] = lowBound[m];
		}
	}
}

static int the_worst_one (Matrix_t* M, int *PF, double* exchv, double* distance, double* angle) {
	int 		rowDim = M->rowDim;
	int 		colDim = M->colDim;
	Matrix_t*	xValue = Matrix_min (M);
	int		isExtreme[rowDim+10];
	int		find[colDim+10];
	int		pf=0, one=-1; 
	double		ex=0, di=0, an=0, d;
	int		i, j;

	memset (isExtreme, 0, rowDim*sizeof(int));
	memset (find, 0, colDim*sizeof(int));
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) if (!find[j]){
			d = M->elements[i*colDim+j] - xValue->elements[j];
			if (d < DBL_EPSILON) {
				isExtreme[i] = 1;
				find[j] = 1;
				break;
			}
		}
	}
	Matrix_free (&xValue);

	for (i=0; i<rowDim; i++) if (!isExtreme[i]) {
		if (PF[i] > pf) {
			one = i;
			pf = PF[i];
			ex = exchv[i];
			di = distance[i];
			an = angle[i];
		} else if (PF[i] == pf && pf == 1) {
			if (exchv[i]*distance[i] < ex*di) {
				one = i;
				pf = PF[i];
				ex = exchv[i];
				di = distance[i];
				an = angle[i];
			}
		} else if (PF[i] == pf && pf > 1) {
			if (angle[i] < an) {
				one = i;
				pf = PF[i];
				ex = exchv[i];
				di = distance[i];
				an = angle[i];
			}
		}
	}
	return one;
}
void Population_remove_one (Population_t *pop) {			// remove one individual by hv  
	int 		i, j, k;
	int 		popSize = pop->var->rowDim;
	int 		numVar  = pop->var->colDim;
	int 		numObj  = pop->obj->colDim;
	double		exchv[popSize+10];
	double		distance[popSize+10];
	double		angle[popSize+10];
	int		PF[popSize+10];
	double		y[numObj+10], t;
	Matrix_t*	Obj  = Matrix_norm (pop->obj); 
	Matrix_t*	M    = NULL;
	Matrix_t*	F    = NULL;
	List_t*		list = ndSort (Obj);
	Link_t*		link = NULL;
	int		*array = NULL;

	// 1. compute exchv
	for (i=0; i<popSize; i++) {
		M = Matrix_dup (Obj);	
		memcpy (y, M->elements+i*numObj, numObj*sizeof(double));
		for (j=0; j<numObj; j++) {
			M->elements[i*numObj+j] = M->elements[(popSize-1)*numObj+j];
		}
		M->rowDim -= 1;
		for (k=0; k<M->rowDim; k++) {
			for (j=0; j<numObj; j++) {
				if (M->elements[k*numObj+j] < y[j]) {
					M->elements[k*numObj+j] = y[j];
				}
			}
		}
		F = Matrix_front (M);
		for (j=0, exchv[i]=1; j<numObj; j++) {
			exchv[i] *= (1.1 - y[j]);
		}
		exchv[i] -= hv (F);
		Matrix_free (&F);
		Matrix_free (&M);
	}
	Matrix_free (&Obj);

	// 2. set the PF
	i=1;
	link=list->list_head;
	while (link) {
		array = Link2Array (link);
		for (j=array[0]; j>0; j--) {
			PF[array[j]] = i;
		}
		free (array);
		i++;
		link=link->next;
	}
	List_free (&list);

	for (i=0; i<popSize; i++) {
		distance[i] = 1.0e+100;
		angle[i] = 1.0e+100;
	}

	// 3. compute distance
	for (i=0; i<popSize; i++) if (1 == PF[i]){
		for (j=i+1; j<popSize; j++) if (1 == PF[j]) {
			t = distance_p2p (pop->obj->elements+i*numObj, pop->obj->elements+j*numObj, numObj); 	
			if (t < distance[i]) {
				distance[i] = t;
			}
			if (t < distance[j]) {
				distance[j] = t;
			}
		}
	}

	// 4. compute angle 
	for (i=0; i<popSize; i++) {
		for (j=i+1; j<popSize; j++) {
			t = vector_angle (pop->obj->elements+i*numObj, pop->obj->elements+j*numObj, numObj); 	
			if (PF[j] <= PF[i] && t < angle[i]) {
				angle[i] = t;
			}
			if (PF[i] <= PF[j] && t < angle[j]) {
				angle[j] = t;
			}
		}
	}
	
	// 5. remove the worst one
	i = the_worst_one (pop->obj, PF, exchv, distance, angle);
	for (j=0; j<numVar; j++) {
		pop->var->elements[i*numVar+j] = pop->var->elements[(popSize-1)*numVar+j];
	}
	for (j=0; j<numObj; j++) {
		pop->obj->elements[i*numObj+j] = pop->obj->elements[(popSize-1)*numObj+j];
	}
	pop->var->rowDim -=1;	
	pop->obj->rowDim -=1;	
}

void Population_optimize_next (Population_t *pop, Optimize_t* opt) {
	int 		i, j, k;
	int		popSize= pop->var->rowDim;
	int		numVar = pop->var->colDim;
	int		numObj = pop->obj->colDim;
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double 		p1[numVar], p2[numVar], c1[numVar], c2[numVar], lb[numVar], ub[numVar];
	int		r1, r2, r3;
	int 		queue[numVar], n=0;	
	int*		array = NULL;
	int*		rank = NULL;
	Matrix_t*	Obj = NULL;
	Matrix_t*	Var = NULL;
	double 		d;

	// get variables
	if (opt->var) {
		array = Link2Array(opt->var);
		memcpy (queue, array+1, array[0]*sizeof (int));
		n = array[0];
		free (array);
	} else {
		for (i=0; i<numVar; i++) {
			queue[i] = i;
		}
		n = numVar;
	}

	// select parents r1, r2
	if (opt->select == Random) {
		r1 = rand() % popSize;
		r2 = rand() % popSize;
		r3 = popSize;
	} else  {
		i = rand() % popSize;
		j = rand() % popSize;
		r1 = (i < j) ? i : j;			// r1
		i = rand() % popSize;
		j = rand() % popSize;
		r2 = (i < j) ? i : j;			// r2
		r3 = popSize;
	}

	// get p1 p2, lb, ub
	for (k=0; k<n; k++) {
		j = queue[k];
		p1[k] = pop->var->elements[r1*numVar+j];
		p2[k] = pop->var->elements[r2*numVar+j];
		lb[k] = lowBound[j];
		ub[k] = uppBound[j];
	}

	// crossover and mutation
	realbinarycrossover(p1, p2, c1, c2, 1.0, n, lb, ub);
       	realmutation(c1, 1.0/n, n, lb, ub);

	// repreduce = Half
	if (Half == opt->reproduce && r1 != r2) for (k=0; k<n; k++) {
		c1[k] = 0.5*(p1[k] + p2[k]);
	}
		
	// construct child with index of r3
	memcpy (pop->var->elements+r3*numVar, pop->var->elements+r1*numVar, numVar*sizeof(double));
	for (k=0; k<n; k++) {
		j=queue[k];
		pop->var->elements[r3*numVar+j] = c1[k];
	}

	// check if the new values of variables are close to repeller
	if (opt->repeller) for (i=0; i<opt->repeller->rowDim; i++) {
		for (j=0; j<numVar; j++) if (pop->I[j]) {
			d = 1.0e-3*(uppBound[j] - lowBound[j]);
			if (fabs(pop->var->elements[r3*numVar+j] - opt->repeller->elements[i*numVar+j]) < d) {
				pop->var->elements[r3*numVar+j] =  lowBound[j] + randu()*(uppBound[j]-lowBound[j]);
			}
		}
	}
	
	// check if execute difference
	if (opt->isExec_difference) {
		Population_exec_difference (pop, pop->var->elements+r3*numVar);
	}

	// evaluate 
	Problem_evaluate (pop->var->elements+r3*numVar,numVar,pop->obj->elements+r3*numObj,numObj);
	
	//
	pop->var->rowDim += 1;
	pop->obj->rowDim += 1;

	// rank
	switch (opt->rank) {
		case Crowding: 
			rank = Matrix_rank_by_crowdingdistance (pop->obj); break;
		case Hypervolume:
			rank = Matrix_rank_by_hypervolume (pop->obj); break;
		case Sum:
			rank = Matrix_rank_by_sum (pop->obj); break;
		default:
			rank = Matrix_rank_by_crowdingdistance (pop->obj);
	}

	Var = Matrix_dup (pop->var);
	Obj = Matrix_dup (pop->obj);
	for (i=0; i<popSize; i++) {
		memcpy (pop->var->elements+i*numVar, Var->elements+rank[i]*numVar, numVar*sizeof (double));
		memcpy (pop->obj->elements+i*numObj, Obj->elements+rank[i]*numObj, numObj*sizeof (double));
	}
	Matrix_free (&Var);
	Matrix_free (&Obj);
	free (rank);

	//
	pop->var->rowDim -= 1;
	pop->obj->rowDim -= 1;
}

Population_t *Population_reproduce (Population_t *P, char type) {
	int i;
	double *p1=NULL, *p2=NULL;
	double *c=NULL, *f=NULL;
	Population_t *Q = Population_dup (P);
	
	for (i=0; i<P->var->rowDim; i++) {
		if (type == 'T') {
			p1 = &P->var->elements[tournament(P->var->rowDim)*P->var->colDim];
			p2 = &P->var->elements[tournament(P->var->rowDim)*P->var->colDim];
		} else {
			p1 = &P->var->elements[(int)(randu()*P->var->rowDim)*P->var->colDim];
			p2 = &P->var->elements[(int)(randu()*P->var->rowDim)*P->var->colDim];
		}

		c = &Q->var->elements[i*P->var->colDim];
		f = &Q->obj->elements[i*P->obj->colDim];
		SBX_reproduce (p1, p2, c, f);
	}
	return Q;
}

Population_t *Population_mutation (Population_t *P) {
	int i;
	double *c=NULL, *f=NULL;
	Population_t *Q = Population_dup (P);
	
	for (i=0; i<P->var->rowDim; i++) {
		c = &Q->var->elements[i*P->var->colDim];
		f = &Q->obj->elements[i*P->obj->colDim];
		SBX_mutation (c, f);
	}
	return Q;
}
Population_t *Population_eliminate (Population_t *pop) {
	Population_t *dup = Population_dup (pop);

	dup->var->rowDim /=2;
	dup->obj->rowDim /=2;

	return dup;
}

Population_t *Population_eliminate (Population_t *pop, int threshold) {
	Population_t *dup = Population_dup (pop);

	dup->var->rowDim = threshold;
	dup->obj->rowDim = threshold;

	return dup;
}

int Population_getSize (Problem_t *problem) {
	int nPop=100;
	Matrix_t *Z = Population_reference (problem);

	nPop = Z->rowDim;
	Matrix_free (&Z);

	return nPop;
} 

int Population_getGen (Problem_t *problem) {

	return 1000;
} 

Matrix_t* Population_reference (Problem_t *problem) {
	Matrix_t *Z=NULL, *T=NULL;

	switch (problem->numObj) {
		case 2: 				
			Z = H1_sample (2, 99);	break;	// nPop = 100 
		case 3: 				
			Z = H1_sample (3, 13);	break;	// nPop = 105   // Z = H1_sample (3, 12);	break;	// nPop = 91
		case 5:
			Z = H1_sample (5, 6);	break; 	// nPop = 210
		case 8:
			Z = H1_sample (8, 3);	
			T = H1_sample (8, 2);	
			Matrix_cat (&Z, T);
			Matrix_free (&T);
			break;			 	// nPop = 156
		case 10:
			Z = H1_sample (10, 3);	
			T = H1_sample (10, 2);	
			Matrix_cat (&Z, T);
			Matrix_free (&T);
			break;				// nPop = 275
			
		case 15:
			Z = H1_sample (15, 2);	
			T = H1_sample (15, 1);	
			Matrix_cat (&Z, T);
			Matrix_free (&T);
			break;				// nPop = 135
		default : 
			fprintf (stderr, "%s:%d: Cannot generate reference points for the problem\n",__FILE__, __LINE__);
			exit (-1);
	}
	return Z;
}

Population_t *Population_Front1 (Population_t *pop) {
	List_t *list = ndSort (pop->obj, 1);
	int *arr = Link2Array (list->list_head);
	Population_t *front = (Population_t *)malloc (sizeof (Population_t));
	
	front->obj = Matrix_sub (pop->obj, arr+1, arr[0]);
	front->var = Matrix_sub (pop->var, arr+1, arr[0]);

	free (arr);
	List_free (&list);

	return front;
}

Population_t *Population_compress (Population_t *pop) {
	Population_t *C = (Population_t *)malloc (sizeof (Population_t));
	
	C->obj = Matrix_compress (pop->obj);
	C->var = Matrix_compress (pop->var);

	return C;
}

