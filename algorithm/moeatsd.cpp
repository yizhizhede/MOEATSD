#include "shape.h"
#include "population.h"
#include "algebra.h"
#include "recombination.h"
#include "terminal.h"
#include "snapshot.h"
#include "myrandom.h"
#include "parameter.h"
#include "dominate.h"
#include "problem.h"
#include "crowding.h"
#include "rank.h"
#include "mystring.h"
#include "interactive.h"
#include "moeatsd.h"
#include "model.h"
#include "bpmlp.h"
#include <string.h>
#include <float.h>

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MAX_NUM_SOLUTIONS 10000
#define OPT_PRINT 0

/*************************************************************************************************************/
/**************  Main framework  *****************************************************************************/
/*************************************************************************************************************/

// Aggregation Function 
static double g_te (double *f, double *lambda, int len);

// SLPSO
static void SLPSO_optimize (Population_t* pop, double (*g_func)(double *f, double *w, int M), double *gW); //SLPSO

// TSD 
static void TSD_analysis_correlation (Population_t* pop);		// correlation between variables and objectives
static void TSD_analysis_landscape (Population_t* pop);
static void TSD_analysis_conflict (Population_t* pop);
static void TSD_conflict_space (Population_t* pop);
static void TSD_satisfy_condition_LS (Population_t*pop, double *x);	// Conditions of landscape
static void TSD_satisfy_condition_CS (Population_t*pop, double *parent, double* child);	//  Conditions of conflict space 
static void TSD_optimize (Population_t* pop);

// Main framewok of MOEA/TSD
Population_t* moeatsd (Problem_t *problem) {
	Population_t* 	pop = NULL;
	int		numVar; 
	int		numObj;
	double		gW[100];

	// 0.1 new a population
	pop 	= Population_new (problem, (Parameter_get())->popSize);
	numVar  = pop->var->colDim;
	numObj  = pop->obj->colDim;

	// 0.2 reallocate memory 
	pop->var->elements = (double *)realloc (pop->var->elements, (MAX_NUM_SOLUTIONS+10)*numVar*sizeof (double));
	pop->obj->elements = (double *)realloc (pop->obj->elements, (MAX_NUM_SOLUTIONS+10)*numObj*sizeof (double));

	// 0.3 alloc memory Conflict variable and Coflict space 
	pop->CV = (int *)calloc (numVar, sizeof (int));
	pop->CS = (int *)calloc (numVar, sizeof (int));

	// 0.4 set gW
	for (int i=0; i<numObj; i++) {gW[i] = 0.01 + randu();}

	// 0.5 print init pop
	isTerminal (pop);

	// 1 correlation 
	TSD_analysis_correlation (pop);
	
	// 2 analysis landscape 
	TSD_analysis_landscape (pop);
	
	// 3 get an optimized soluton
	SLPSO_optimize (pop, g_te, gW); 

	// 4 analysis conflict variable
	TSD_analysis_conflict (pop);

	// print CV
	printf ("CV = { "); for (int i=0; i<numVar; i++) if (1 == pop->CV[i]) printf ("%d ", i); printf ("}\n");

	// 5 set conflict space
	TSD_conflict_space (pop);

	// print CS
	printf ("CS = { "); for (int i=0; i<numVar; i++) if (1 == pop->CS[i]) printf ("%d ", i); printf ("}\n");

	// 6 TSD optimize
	TSD_optimize (pop);

	// return
	return pop;
}

/*************************************************************************************************************/
/*************** TSD optimize ********************************************************************************/
/*************************************************************************************************************/

typedef struct TSD_Pop_tag {
	double*	var;
	double* vel;
	double* obj;
	int	popSize;
	int 	numVar;
	int	numObj;	
} TSD_Pop_t;

static TSD_Pop_t* TSD_Pop_new (int popSize, int numVar, int numObj) {
	TSD_Pop_t* pop = NULL;

	pop = (TSD_Pop_t *)calloc (1, sizeof (TSD_Pop_t));
	pop->var = (double *)calloc (popSize*numVar+10, sizeof (double));
	pop->vel = (double *)calloc (popSize*numVar+10, sizeof (double));
	pop->obj = (double *)calloc (popSize*numObj+10, sizeof (double));
	pop->popSize = popSize;
	pop->numObj  = numObj;
	pop->numVar  = numVar;
	return pop;
}

static void TSD_Pop_free (TSD_Pop_t* pop) {
	free (pop->obj);
	free (pop->vel);
	free (pop->var);
	free (pop);
}

static void TSD_init_pop (Population_t* pop, TSD_Pop_t* popOffspring) {
	double*	lowBound = Problem_getLowerBound ();
	double*	uppBound = Problem_getUpperBound ();
	int 	popSize  = popOffspring->popSize;
	int	numVar   = popOffspring->numVar;
	int	numObj 	 = popOffspring->numObj;
	int 	i, j;

	// var and obj
	memcpy (popOffspring->var, pop->var->elements, pop->var->rowDim*numVar*sizeof (double));
	memcpy (popOffspring->obj, pop->obj->elements, pop->obj->rowDim*numObj*sizeof (double));
	for (i=pop->var->rowDim; i<popSize; i++) {
		for (j=0; j<numVar; j++) {
			popOffspring->var[i*numVar+j] = lowBound[j] + randu()*(uppBound[j]-lowBound[j]);
		}
		Problem_evaluate (popOffspring->var+i*numVar, numVar, popOffspring->obj+i*numObj, numObj);
	}

	// velocity
	memset (popOffspring->vel, 0, popSize*numVar*sizeof (double));
}

static double getTheta_fitness;
static double TSD_get_theta () {
	double 	t = 0;

	if (!getTheta_fitness) { 
		getTheta_fitness = (double)Problem_getFitness ();
	}
	t = (Problem_getFitness () - getTheta_fitness)/(Problem_getLifetime () - getTheta_fitness);

	return  t * t;
}

static void TSD_get_gamma (Matrix_t *M, double *gamma) {
	int 	rowDim = M->rowDim;
	int	colDim = M->colDim;
	double 	distance = 0;
	int 	i, j;

	for (i=0; i<rowDim; gamma[i] = 1.0e+100, i++) {};
	for (i=0; i<rowDim-1; i++) {
		for (j=i+1; j<rowDim; j++) {
			distance = distance_p2p (M->elements+i*colDim, M->elements+j*colDim, colDim);
			if (distance < gamma[i]) 
				gamma[i] = distance;
			if (distance < gamma[j]) 
				gamma[j] = distance;
		}
	}
}

static void TSD_get_associate (Matrix_t* M, Matrix_t* V, int *associate) {
	int 	i, j, b;
	double	minDis, dis;

	// associate
	for (i=0; i<M->rowDim; i++) {
		b = 0;
		minDis = distance_p2p (M->elements+i*M->colDim, V->elements+0*V->colDim, V->colDim);	
		for (j=1; j<V->rowDim; j++) {
			dis = distance_p2p (M->elements+i*M->colDim, V->elements+j*V->colDim, V->colDim);	
			if (dis < minDis) {
				minDis = dis;
				b = j;
			}
		}
		associate[i] = b;
	}
}

static Matrix_t* SBV_get_Var (TSD_Pop_t* popOffspring, int* byVar, int len) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		popSize	 = popOffspring->popSize;
	int 		numVar   = popOffspring->numVar;
	Matrix_t*	Var = Matrix_new (popSize, len);
	int		i, j, k;

	// get Var
	for (k=0; k<len; k++) {
		j = byVar[k];
		for (i=0; i<popSize; i++) {
			Var->elements[i*len+k] = (popOffspring->var[i*numVar+j]-lowBound[j])/(uppBound[j]-lowBound[j]);
		}
	}
	return Var;
}

static void TSD_selection_by_var (TSD_Pop_t* popOffspring, Matrix_t *V, double theta, double* gammaV, int *byVar, int *sel, TSD_Pop_t* popSelectbyVar) {
	int		popSize	= popOffspring->popSize;
	int 		numVar  = popOffspring->numVar;
	int 		numObj 	= popOffspring->numObj;
	Matrix_t*	Var 	= NULL;
	int 		associateV[popSize+10];
	int 		vis[V->rowDim+10];
	int		i, j, a, b, n=0;
	double		DPD, distance, length, minDPD;

	// get Var
	Var = SBV_get_Var (popOffspring, byVar, V->colDim);
	
	// associate
	TSD_get_associate (Var, V, associateV);

	// set vis
	memset (vis, 0, V->rowDim*sizeof (int));

	// select one for each reference point
	for (a=0, n=0; a<popSize; a++) {
		i = associateV[a];	// reference index
		if (0 == vis[i]) {
			vis[i] = 1;
			b = -1;
			for (j=0; j<popSize; j++) if (associateV[j] == i) { // individual index
				distance = distance_p2p (Var->elements+j*Var->colDim, V->elements+i*V->colDim, V->colDim);
				length = norm (popOffspring->obj+j*numObj, numObj);
				DPD = (1 + numObj * theta * distance / gammaV[i]) * length;
				if (-1 == b) {
					minDPD = DPD;
					b = j;
				} else if (DPD < minDPD) {
					minDPD = DPD;
					b = j;
				}
			}
			if (b != -1 && 0 == sel[b]) {
				memcpy (popSelectbyVar->var+n*numVar, popOffspring->var+b*numVar, numVar*sizeof (double));
				memcpy (popSelectbyVar->vel+n*numVar, popOffspring->vel+b*numVar, numVar*sizeof (double));
				memcpy (popSelectbyVar->obj+n*numObj, popOffspring->obj+b*numObj, numObj*sizeof (double));
				sel[b] = 1;	// set flag of selection
				n++;	
			}
		}
	}
	
	// Set size
	popSelectbyVar->popSize = n;

	// free Var
	Matrix_free (&Var);
}

//
static Matrix_t* SBO_get_obj (TSD_Pop_t* popOffspring) {
	int 		popSize	= popOffspring->popSize;
	int		numObj  = popOffspring->numObj;
	Matrix_t*	Obj   	= Matrix_new (popSize, numObj);
	double		minObj, maxObj, L;
	int		i, j;

	for (j=0; j<numObj; j++) {
		minObj = popOffspring->obj[0*numObj+j]; 
		maxObj = popOffspring->obj[0*numObj+j]; 
		for (i=1; i<popSize; i++) {
			if (popOffspring->obj[i*numObj+j] < minObj) 
				minObj = popOffspring->obj[i*numObj+j];
			if (popOffspring->obj[i*numObj+j] > maxObj) 
				maxObj = popOffspring->obj[i*numObj+j];
		}
		L = maxObj - minObj;
		if (L > DBL_EPSILON) for (i=0; i<popSize; i++) {
			Obj->elements[i*numObj+j] = (popOffspring->obj[i*numObj+j] - minObj) / L;
		} else for (i=0; i<popSize; i++) {
			Obj->elements[i*numObj+j] = popOffspring->obj[i*numObj+j] - minObj;
		}
	}

	//
	for (i=0; i<popSize; i++) {
		L = sum (Obj->elements+i*numObj, numObj);	
		if (L > DBL_EPSILON) for (j=0; j<numObj; j++) {
			Obj->elements[i*numObj+j] /= L;
		}
	}
	return Obj;
}

static void TSD_selection_by_obj (TSD_Pop_t* popOffspring, Matrix_t* Z, double theta, double* gammaZ, int* sel, TSD_Pop_t* popSelectbyObj) {
	int 		popSize	= popOffspring->popSize;
	int		numVar  = popOffspring->numVar;
	int		numObj  = popOffspring->numObj;
	Matrix_t*	Obj 	= NULL;
	int		associateZ[popSize+10];
	int 		vis[Z->rowDim+10];
	int 		i, j, a, b, n = 0;
	double		DPD, distance, length, minDPD;
	
	// get Obj
	Obj = SBO_get_obj (popOffspring);

	// associate
	TSD_get_associate (Obj, Z, associateZ);

	// set vis 
	memset (vis, 0, Z->rowDim*sizeof (int));

	// select one for each reference point
	for (a=0, n=0; a<popSize; a++) {
		i = associateZ[a];	// reference vector index
		if (0 == vis[i]) {
			vis[i] = 1;
			b = -1;
			for (j=0; j<popSize; j++) if (associateZ[j] == i) { 	// solution index
				distance = distance_p2p (Obj->elements+j*numObj, Z->elements+i*numObj, numObj);
				length = norm (popOffspring->obj+j*numObj, numObj);
				DPD = (1 + numObj * theta * distance / gammaZ[i]) * length;
				if (-1 == b) {
					minDPD = DPD;
					b = j;
				} else if (DPD < minDPD) {
					minDPD = DPD;
					b = j;		
				}
			}
			if (b != -1 && 0 == sel[b]) {
				memcpy (popSelectbyObj->var+n*numVar, popOffspring->var+b*numVar, numVar*sizeof (double));
				memcpy (popSelectbyObj->vel+n*numVar, popOffspring->vel+b*numVar, numVar*sizeof (double));
				memcpy (popSelectbyObj->obj+n*numObj, popOffspring->obj+b*numObj, numObj*sizeof (double));
				sel[b] = 1;	// set flag of selection
				n++;	
			}
		}
	}

	// Set size
	popSelectbyObj->popSize = n;

	// free Obj
	Matrix_free (&Obj);
}

//
typedef struct TSD_EP02_tag {
	double*	var;	
	double*	obj;	
	int*	used;		// a array
	int	len;
	int*	idle;		// a stack
	int	top;
	int*	uppNeighbor;
	int*	lowNeighbor;
} TSD_EP02_t;

//
static TSD_EP02_t* TSD_EP02_new (int EPSize, int numVar, int numObj) {
	TSD_EP02_t*	tsdEP 	= NULL;
	int		i, j;

	// allocation	
	EPSize += 10;
	tsdEP = (TSD_EP02_t*)calloc (1, sizeof (TSD_EP02_t));
	tsdEP->var = (double*) calloc (EPSize*numVar+10, sizeof (double));
	tsdEP->obj = (double*) calloc (EPSize*numObj+10, sizeof (double));
	tsdEP->used = (int*) calloc (EPSize+10, sizeof (int));
	tsdEP->idle = (int*) calloc (EPSize+10, sizeof (int));
	tsdEP->uppNeighbor = (int*) calloc (EPSize*numObj+10, sizeof (int));
	tsdEP->lowNeighbor = (int*) calloc (EPSize*numObj+10, sizeof (int));

	// initialization
	tsdEP->len = 0;
	tsdEP->top = 0;
	for (i=2; i<EPSize; tsdEP->idle[tsdEP->top++]=i, i++) {};
	for (j=0; j<numObj; j++) {
		tsdEP->obj[0*numObj+j] = -1.0e+100;
		tsdEP->obj[1*numObj+j] =  1.0e+100;
		tsdEP->uppNeighbor[0*numObj+j] = 1;
		tsdEP->lowNeighbor[0*numObj+j] = 0;
		tsdEP->uppNeighbor[1*numObj+j] = 1;
		tsdEP->lowNeighbor[1*numObj+j] = 0;
	}
	return tsdEP;
}

//
static void TSD_EP02_add (Population_t* pop, TSD_EP02_t* tsdEP, double *var, double* obj) {
	int	numVar	 = pop->var->colDim;
	int	numObj   = pop->obj->colDim;
	int	i, j, k, low;

	// pop one from idle
	i = tsdEP->idle[tsdEP->top-1];
	tsdEP->top--;

	// put one into used
	tsdEP->used[tsdEP->len] = i;
	tsdEP->len++;

	// copy to i
	memcpy (tsdEP->var+i*numVar, var, numVar*sizeof (double));
	memcpy (tsdEP->obj+i*numObj, obj, numObj*sizeof (double));

	// insert into link 
	for (j=0; j<numObj; j++) {
		k = tsdEP->uppNeighbor[0*numObj+j];
		while (tsdEP->obj[k*numObj+j] < obj[j]) {
			k = tsdEP->uppNeighbor[k*numObj+j];
		}
		low = tsdEP->lowNeighbor[k*numObj+j];
		tsdEP->uppNeighbor[low*numObj+j] = i;
		tsdEP->lowNeighbor[i*numObj+j] = low;
		tsdEP->uppNeighbor[i*numObj+j] = k;
		tsdEP->lowNeighbor[k*numObj+j] = i;
	}
}

//
static void TSD_EP02_delete (Population_t* pop, TSD_EP02_t* tsdEP, int i) {
	int	numObj   = pop->obj->colDim;
	int	upp, low;
	int	j, a;

	// push i into idle
	tsdEP->idle[tsdEP->top] = i;
	tsdEP->top++;

	// delete i from used
	for (a=0; a<tsdEP->len; a++) if (i == tsdEP->used[a]) {
		tsdEP->used[a] = tsdEP->used[tsdEP->len-1];
		tsdEP->len--;
		break;
	}

	// delete link 
	for (j=0; j<numObj; j++) {
		upp = tsdEP->uppNeighbor[i*numObj+j];
		low = tsdEP->lowNeighbor[i*numObj+j];
		tsdEP->lowNeighbor[upp*numObj+j] = low;
		tsdEP->uppNeighbor[low*numObj+j] = upp;
	}
}

//
#define TUEP_LEN 200
static double TUEP_historyObj[1000*100];
static int    TUEP_cur;
//
static int TUEP_update_history (double* obj, int numObj) {
	int	i, j;
	double	t;

	for (i=TUEP_cur-1; i>=0; i--) {
		for (j=0; j<numObj; j++) {
			t = TUEP_historyObj[i*numObj+j] - obj[j];
			if (t > DBL_EPSILON || t < -DBL_EPSILON) {
				break;
			}
		}
		if (j >= numObj) {
			return 0;
		}
	}
	for (i=TUEP_LEN-1; i>=TUEP_cur; i--) {
		for (j=0; j<numObj; j++) {
			t = TUEP_historyObj[i*numObj+j] - obj[j];
			if (t > DBL_EPSILON || t < -DBL_EPSILON) {
				break;
			}
		}
		if (j >= numObj) {
			return 0;
		}
	}
	memcpy (TUEP_historyObj+TUEP_cur*numObj, obj, numObj*sizeof (double));
	TUEP_cur = (TUEP_cur + 1) % TUEP_LEN;
	return 1;
}

static int TEU_less (double *A, double *B, int len) {
	int 	i;

	for (i=0; i<len; i++) {
		if (A[i] < B[i]) { 
			return 1;
		} else if (A[i] == B[i]) {
			continue;
		} else {
			return 0;
		}
	}
	return 0;
}

//
static void TSD_EP02_update (Population_t* pop, TSD_EP02_t* tsdEP, double *var, double* obj) {
	int	numObj   = pop->obj->colDim;
	double	crowdingD[Parameter_get()->popSize+20];
	double	minCD = 0, minV, maxV;
	double	minObj[numObj+10];
	double	buff[numObj+10];
	int 	dominate = 0; 
	int	upp, low;
	int 	i, j, k, a, b;

	// update history
	if (0 == TUEP_update_history (obj, numObj)) return;

	// compare one with EP
	for (a=0; a<tsdEP->len; a++) {
		i = tsdEP->used[a];
		dominate = isDominate (tsdEP->obj+i*numObj, obj, numObj);
		if (1 == dominate || 0 == dominate) {
			return;
		}
	}

	// delete one dominated by i
	for (a=0; a<tsdEP->len; a++) {
		i = tsdEP->used[a];
		dominate = isDominate (obj, tsdEP->obj+i*numObj, numObj);
		if (1 == dominate) {
			TSD_EP02_delete (pop, tsdEP, i);// delete one from EP
			a--;				// check 'a' again
		}
	}

	// add one to EP
	TSD_EP02_add (pop, tsdEP, var, obj);	

	// if size is smaller
	if (tsdEP->len <= Parameter_get()->popSize) return;

	// init crowdingD
	memset (crowdingD, 0, (Parameter_get()->popSize+20)*sizeof (double));

	// conpute crowdinD
	for (i=0; i<numObj; i++) {
		for (j=0; j<numObj; j++) {
			minObj[(j+i) % numObj] = 1.0e+100;
			b = 1;
		}
		for (a=0; a<tsdEP->len; a++) {
			k = tsdEP->used[a];
			for (j=0; j<numObj; j++) {
				buff[(j+i) % numObj] = tsdEP->obj[k*numObj+j];
			}
			if (1 == TEU_less (buff, minObj, numObj)) {
				memcpy (minObj, buff, numObj*sizeof (double));
				b = k;
			}
		}
		crowdingD[b] = 1.0e+100;
	}
	
	// conpute crowdinD
	for (j=0; j<numObj; j++) {
		maxV = tsdEP->obj[tsdEP->lowNeighbor[1*numObj+j]*numObj+j];
		minV = tsdEP->obj[tsdEP->uppNeighbor[0*numObj+j]*numObj+j];
		if (maxV - minV < DBL_EPSILON) continue;
		i = tsdEP->uppNeighbor[0*numObj+j];
		i = tsdEP->uppNeighbor[i*numObj+j];
		while (1 != tsdEP->uppNeighbor[i*numObj+j]) {
			upp = tsdEP->uppNeighbor[i*numObj+j];
			low = tsdEP->lowNeighbor[i*numObj+j];
			crowdingD[i] += (tsdEP->obj[upp*numObj+j] - tsdEP->obj[low*numObj+j]) / (maxV - minV); 
			i = tsdEP->uppNeighbor[i*numObj+j];
		}
	}

	// find one with minCD
	for (a=0, b=-1; a<tsdEP->len; a++) {
		i = tsdEP->used[a];
		if (-1 == b) {
			minCD = crowdingD[i];
			b = i;
		} else if (crowdingD[i] < minCD) {
			minCD = crowdingD[i];
			b = i;
		}
	}
	
	// delete one from EP
	if (-1 != b) {
		TSD_EP02_delete (pop, tsdEP, b);	
	}
}

//
static void TSD_EP02_update (Population_t* pop, TSD_EP02_t* tsdEP, TSD_Pop_t* popSelectbyVar, TSD_Pop_t* popSelectbyObj) {
	int 	popSize1 = popSelectbyVar->popSize;
	int 	popSize2 = popSelectbyObj->popSize;
	int	numVar	 = pop->var->colDim;
	int	numObj   = pop->obj->colDim;
	int	i, a;

	// update EP using popSelectbyVar
	for (i=0; i<popSize1; i++) {
		TSD_EP02_update (pop, tsdEP, popSelectbyVar->var+i*numVar, popSelectbyVar->obj+i*numObj);
	}

	// update EP using popSelectbyObj
	for (i=0; i<popSize2; i++) {
		TSD_EP02_update (pop, tsdEP, popSelectbyObj->var+i*numVar, popSelectbyObj->obj+i*numObj);
	}

	// copy EP to pop
	for (a=0; a<tsdEP->len; a++) {
		i = tsdEP->used[a];
		memcpy (pop->var->elements+a*numVar, tsdEP->var+i*numVar, numVar*sizeof (double));
		memcpy (pop->obj->elements+a*numObj, tsdEP->obj+i*numObj, numObj*sizeof (double));
	}
	
	// set size of pop
	pop->var->rowDim = tsdEP->len;	
	pop->obj->rowDim = tsdEP->len;	
}

typedef struct TSD_EP03_tag {
	Matrix_t*	Z;	//	
	double*		gammaZ;	//		
	double*		minObj;	//
	double*		L;	//
} TSD_EP03_t;

static TSD_EP03_t* TSD_EP03_new (Matrix_t* Z, double *gammaZ) {
	TSD_EP03_t*	tsdEP = NULL;
	int		colDim = Z->colDim;
	int		i;

	// allocate
	tsdEP = (TSD_EP03_t*)calloc (1, sizeof (TSD_EP03_t));
	tsdEP->minObj = (double*)calloc (colDim+10, sizeof (double));
	tsdEP->L      = (double*)calloc (colDim+10, sizeof (double));

	// init
	tsdEP->Z = Z;
	tsdEP->gammaZ = gammaZ;
	for (i=0; i<colDim; i++) {tsdEP->L[i] = 1.0;}

	return tsdEP;
}

static Matrix_t* EP03_get_obj (TSD_Pop_t* popAll, TSD_EP03_t* tsdEP) {
	int 		popSize	= popAll->popSize;
	int		numObj  = popAll->numObj;
	Matrix_t*	Obj   	= Matrix_new (popSize, numObj);
	int		i, j;
	double		L;

	// update minObj
	for (j=0; j<numObj; j++) {
		tsdEP->minObj[j] = popAll->obj[j];
		for (i=1; i<popSize; i++) if (popAll->obj[i*numObj+j] < tsdEP->minObj[j]) {
	 		tsdEP->minObj[j] = popAll->obj[i*numObj+j];
		}
	}

	//
	for (j=0; j<numObj; j++) {
		if (tsdEP->L[j] > DBL_EPSILON) for (i=0; i<popSize; i++) {
			Obj->elements[i*numObj+j] = (popAll->obj[i*numObj+j] - tsdEP->minObj[j]) / tsdEP->L[j];
		} else for (i=0; i<popSize; i++) {
			Obj->elements[i*numObj+j] = popAll->obj[i*numObj+j] - tsdEP->minObj[j];
		}
	}

	//
	for (i=0; i<popSize; i++) {
		L = sum (Obj->elements+i*numObj, numObj);	
		if (L > DBL_EPSILON) for (j=0; j<numObj; j++) {
			Obj->elements[i*numObj+j] /= L;
		}
	}
	return Obj;
}

static int TUL_cmp(const void* a, const void* b) {
	double d = *((double*)a) - *((double*)b);
	if (d > 0) return 1;
	else       return 0;
}
static void TEU_update_L (Population_t* pop, TSD_EP03_t* tsdEP) {
	int	popSize= pop->obj->rowDim;
	int	numObj = pop->obj->colDim;
	double	buff[popSize+10];
	double	L[numObj+10], minL;
	int	k1, k2;
	int	i, j; 

	k1 = popSize / 10;
	k2 = 9*popSize / 10;
	for (j=0; j<numObj; j++) {
		for (i=0; i<popSize; i++) {
			buff[i] = pop->obj->elements[i*numObj+j];
		}
		qsort (buff, popSize, sizeof (double), TUL_cmp);
		L[j] = buff[k2] - buff[k1];
	}

	minL=L[0];
	for (j=1; j<numObj; j++) if (L[j] < minL) {
		minL = L[j];
	}
	for (j=0; j<numObj; j++) {
		tsdEP->L[j] =  (int)(L[j]/minL);
	}
}

//
static void TSD_EP03_update (Population_t* pop, TSD_EP03_t* tsdEP, TSD_Pop_t* popSelectbyVar, TSD_Pop_t* popSelectbyObj) {
	int	 	popSize0 = pop->var->rowDim;
	int 		popSize1 = popSelectbyVar->popSize;
	int 		popSize2 = popSelectbyObj->popSize;
	int		numVar	 = pop->var->colDim;
	int		numObj	 = pop->obj->colDim;
	Matrix_t*	Obj 	 = NULL;
	TSD_Pop_t*	popAll 	 = TSD_Pop_new (popSize0+popSize1+popSize2, numVar, numObj);
	Matrix_t*	Z 	 = tsdEP->Z;
	double*		gammaZ 	 = tsdEP->gammaZ;
	int		associateZ[popAll->popSize+10];
	int		sel[popAll->popSize+10];
	int		vis[Z->rowDim+10];
	int		i, j, k, a, b, n = 0;
	double		distance, length, DPD, minDPD; 
	double		SDE[popAll->popSize+10], maxSDE, t, d;

	// popSelectbyVar
	for (i=0, n=0; i<popSize1; i++) if (1 == TUEP_update_history (popSelectbyVar->obj+i*numObj, numObj)) {
		memcpy(popAll->var+n*numVar, popSelectbyVar->var+i*numVar, numVar*sizeof (double));
		memcpy(popAll->obj+n*numObj, popSelectbyVar->obj+i*numObj, numObj*sizeof (double));
		n++;
	}
	// popSelectbyVar
	for (i=0; i<popSize2; i++) if (1 == TUEP_update_history (popSelectbyObj->obj+i*numObj, numObj)) {
		memcpy(popAll->var+n*numVar, popSelectbyObj->var+i*numVar, numVar*sizeof (double));
		memcpy(popAll->obj+n*numObj, popSelectbyObj->obj+i*numObj, numObj*sizeof (double));
		n++;
	}
	if (n < 1) {
		TSD_Pop_free (popAll);
		return;
	}

	// copy pop to popAll
	memcpy(popAll->var+n*numVar, pop->var->elements, popSize0*numVar*sizeof (double));
	memcpy(popAll->obj+n*numObj, pop->obj->elements, popSize0*numObj*sizeof (double));

	// set popAll size
	popAll->popSize = n + popSize0;

	// get Obj
	Obj =  EP03_get_obj (popAll, tsdEP);

	// associate
	TSD_get_associate (Obj, Z, associateZ);

	// set vis and sel
	memset (vis, 0, Z->rowDim*sizeof (int));
	memset (sel, 0, popAll->popSize*sizeof (int));

	// select one for each reference point
	for (a=0, n=0; a<popAll->popSize; a++) {
		i = associateZ[a];	// reference vector index
		if (0 == vis[i]) {
			vis[i] = 1;
			b = -1;
			for (j=0; j<popAll->popSize; j++) if (associateZ[j] == i) { 	// solution index
				distance = distance_p2p (Obj->elements+j*numObj, Z->elements+i*numObj, numObj);
				length = norm (popAll->obj+j*numObj, numObj);
				DPD = (1 + numObj * distance / gammaZ[i]) * length;
				if (-1 == b) {
					minDPD = DPD;
					b = j;
				} else if (DPD < minDPD) {
					minDPD = DPD;
					b = j;		
				}
			}
			if (b != -1) {
				sel[b] = 1;	// set flag of selection
				n++;
			}
		}
	}

	// compute SDE
	if (n < Z->rowDim && n < popAll->popSize) {
		for (i=0; i<popAll->popSize; i++) if (0 == sel[i]) {
			SDE[i] = 1.0e+100;
			for (k=0; k<popAll->popSize; k++) if (1 == sel[k]) {
				for (j=0, t=0; j<numObj; j++) {
					d = (popAll->obj[k*numObj+j] - popAll->obj[i*numObj+j]);
					if (d > 0) {
						t += d * d;
					}
				}
				// t = sqrt (t); replace by t = t;
				if (t < SDE[i] ) {
					SDE[i] = t;
				}
			}
		}
	}

	// fill empty 
	while (n < Z->rowDim && n < popAll->popSize) {
		maxSDE = -1.0e+100;
		b = -1;
		for (i=0; i<popAll->popSize; i++) if (0 == sel[i]) {
			if (SDE[i] > maxSDE) {
				maxSDE = SDE[i];
				b = i;
			}
		}
		if (-1 != b) {
			sel[b] = 1;	// set flag of selection
			n++;
			// update SDE
			for (i=0; i<popAll->popSize; i++) if (0 == sel[i]) {
				for (j=0, t=0; j<numObj; j++) {
					d = (popAll->obj[b*numObj+j] - popAll->obj[i*numObj+j]);
					if (d > 0) {
						t += d * d;
					}
				}
				// t = sqrt (t); replace by t = t;
				if (t < SDE[i] ) {
					SDE[i] = t;
				}
			}
		}
	}

	// copy popAll to pop
	for (i=0, k=0; i<popAll->popSize; i++) if (1 == sel[i]) {
		memcpy (pop->var->elements+k*numVar, popAll->var+i*numVar, numVar*sizeof (double));
		memcpy (pop->obj->elements+k*numObj, popAll->obj+i*numObj, numObj*sizeof (double));
		k++;
	}

	// Set size
	pop->var->rowDim = k;
	pop->obj->rowDim = k;

	// update tsdEP
	TEU_update_L (pop, tsdEP);

	// free Obj, popAll
	Matrix_free (&Obj);
	TSD_Pop_free (popAll);
}

static  Matrix_t* TSD_get_RV (int rowDim, int colDim) {
	Matrix_t* 	RV = Matrix_new (rowDim, colDim);
	Matrix_t*	Q  = Matrix_new (2*rowDim, colDim);
	double		Dis[4*rowDim*rowDim+10];
	int		vis[2*rowDim+10];
	int		i, j, k, b, p1, p2;
	int		maxGen = 1000, gen = 0;
	double		child1[colDim+10], child2[colDim+10];
	double		low[colDim+10], upp[colDim+10];
	double		minDis, maxDis, t;

	// 
	if (colDim == 1) {
		for (i=0; i<rowDim; i++) {
			RV->elements[i] = i / (rowDim - 1.0);
		}
		Matrix_free (&Q);
		return RV;
	} 

	// init
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			RV->elements[i*colDim+j] = randu ();
		}
	}

	// set boundary
	for (j=0; j<colDim; j++) {
		low[j] = 0;
		upp[j] = 1;
	}

	// loop
	while (gen < maxGen) {

		// reproduce
		for (i=0; i<rowDim; i++) {
			p1 = rand () % rowDim;
			p2 = rand () % rowDim;
			realbinarycrossover(RV->elements+p1*colDim, RV->elements+p2*colDim, child1, child2, 
				1.0, colDim, low, upp);
			realmutation(child1, 1.0/colDim, colDim, low, upp);
			for (j=0; j<colDim; j++) {
				Q->elements[i*colDim+j] = child1[j];
			}
		}
		memcpy (Q->elements+rowDim*colDim, RV->elements, rowDim*colDim*sizeof (double));

		// distance
		for (i=0; i<2*rowDim; i++) {
			Dis[i*2*rowDim+i] = 0;
			for (j=i+1; j<2*rowDim; j++) {
				Dis[i*2*rowDim+j] = distance_p2p (Q->elements+i*colDim, Q->elements+j*colDim, colDim);
				Dis[j*2*rowDim+i] = Dis[i*2*rowDim+j];
			}
		}
		
		// set vis
		memset (vis, 0, 2*rowDim*sizeof (int));
		p1 = rand () % (2*rowDim);
		vis[p1] = 1;
		
		// select
		for (k=1; k<rowDim ;k++) {
			maxDis = -1.0e+100;
			b = -1;
			for (i=0; i<2*rowDim; i++) if (0 == vis[i]) {
				minDis = 1.0e+100;
				for (j=0; j<2*rowDim; j++) if (1 == vis[j]) {
					t = Dis[i*2*rowDim+j]; 
					if (t < minDis) {
						minDis = t;
					}
				}
				if (minDis > maxDis) {
					maxDis = minDis;
					b = i;
				}
			}
			if (-1 != b) {
				vis[b] = 1;
			}
		}
		
		// 
		for (i=0, k=0; i<2*rowDim; i++) if (1 == vis[i]) {
			memcpy (RV->elements+k*colDim, Q->elements+i*colDim, colDim*sizeof (double));
			k++;
		}
		
		//
		gen++;
	}
	
	Matrix_free (&Q);
	return RV;
}

// reproduce using CSO
static void TSD_reproduce_CSO (Population_t* pop, TSD_Pop_t* popParent, TSD_Pop_t* popOffspring) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int 		popSize	= popParent->popSize;
	int		numVar  = popParent->numVar;
	int		numObj  = popParent->numObj;
	Matrix_t*	M   = Matrix_new (popSize, numObj);
	Matrix_t*	Obj = NULL;
	int*		I   = NULL;
	int 		queue[popSize+10], n=0;
	int		loser, winner;
	double		r1, r2;
	double		angle, minAngle;
	int 		i, j, k = 0, a, b = 0;

	// get Obj
	memcpy (M->elements, popParent->obj, popSize*numObj*sizeof (double));
	Obj = Matrix_norm (M);

	// get rank
	I = rank_by_density (Obj);

	// set queue
	for (i=0, n=0; i<popSize; i++) { 
		queue[n++] = i;
	}

	// loop
	while (n > 0) {
		// get loser and winner
		if (n == 1) {
			loser  = queue[0];
			winner = queue[0]; 
			n--;
		} else {
			// loser
			i = rand() % n;
			loser = queue[i]; 
			queue[i] = queue[n-1];
			n--;
			
			// winner
			for (a=0, minAngle=1.0e+100; a<n; a++) {
				i = queue[a];
				angle = vector_angle (Obj->elements+loser*numObj, Obj->elements+i*numObj, numObj); 
				if (angle < minAngle) {
					minAngle = angle;
					b = a;
				}
			}
			//
			winner = queue[b]; 
			queue[b] = queue[n-1];
			n--;

			// make sure winner is better 
			for (i=0; i<popSize; i++) if (I[i] == winner) break; 
			for (j=0; j<popSize; j++) if (I[j] == loser) break; 
			if (i > j) {
				a = loser; loser = winner; winner = a;	
			}
		}

		// r1, r2
		r1 = randu ();
		r2 = randu ();

		// update loser
		for (j=0; j<numVar; j++) {
			popOffspring->vel[k*numVar+j] = r1*popParent->vel[loser*numVar+j] 
				+ r2*(popParent->var[winner*numVar+j] - popParent->var[loser*numVar+j]);
			popOffspring->var[k*numVar+j] = popParent->var[loser*numVar+j] + popOffspring->vel[k*numVar+j] 
				+ r1*(popOffspring->vel[k*numVar+j] - popParent->vel[loser*numVar+j]);

			//
			if (popOffspring->var[k*numVar+j] > uppBound[j]) {
				popOffspring->var[k*numVar+j] = uppBound[j];
			}
			if (popOffspring->var[k*numVar+j] < lowBound[j]) {
				popOffspring->var[k*numVar+j] = lowBound[j];
			}
		}

		// mutation loser
		realmutation(popOffspring->var+k*numVar, 1.0/numVar, numVar, lowBound, uppBound);

		// satisfy condtion of conflict space
		TSD_satisfy_condition_CS (pop, popParent->var+loser*numVar, popOffspring->var+k*numVar);

		// evaluate loser
		Problem_evaluate (popOffspring->var+k*numVar, numVar, popOffspring->obj+k*numObj, numObj);
		k++;
		
		// copy winner 
		memcpy (popOffspring->vel+k*numVar, popParent->vel+winner*numVar, numVar*sizeof (double));
		memcpy (popOffspring->var+k*numVar, popParent->var+winner*numVar, numVar*sizeof (double));

		// mutation winner
		realmutation(popOffspring->var+k*numVar, 1.0/numVar, numVar, lowBound, uppBound);

		// satisfy condtion of conflict space
		TSD_satisfy_condition_CS (pop, popParent->var+loser*numVar, popOffspring->var+k*numVar);

		// evaluate winner 
		Problem_evaluate (popOffspring->var+k*numVar, numVar, popOffspring->obj+k*numObj, numObj);
		k++;
	}

	// copy popParent to popOffspring
	memcpy (popOffspring->vel+k*numVar, popParent->vel, popSize*numVar*sizeof (double));
	memcpy (popOffspring->var+k*numVar, popParent->var, popSize*numVar*sizeof (double));
	memcpy (popOffspring->obj+k*numObj, popParent->obj, popSize*numObj*sizeof (double));
	
	// set size
	popOffspring->popSize = k + popSize;

	// free
	free (I);
	Matrix_free (&M);
	Matrix_free (&Obj);
}


// TSD reproduce
static void TSD_reproduce (Population_t* pop, TSD_Pop_t* popSelectbyVar, TSD_Pop_t* popSelectbyObj, TSD_Pop_t* popOffspring){
	int 	popSize1= popSelectbyVar->popSize;
	int 	popSize2= popSelectbyObj->popSize;
	int 	popSize	= popSelectbyVar->popSize + popSelectbyObj->popSize;
	int	numVar  = pop->var->colDim;
	int	numObj  = pop->obj->colDim;

	// popParent
	TSD_Pop_t* popParent = TSD_Pop_new (popSize, numVar, numObj);

	// and popSelectbyVar to popParent
	memcpy (popParent->var, popSelectbyVar->var, popSize1*numVar*sizeof (double));
	memcpy (popParent->vel, popSelectbyVar->vel, popSize1*numVar*sizeof (double));
	memcpy (popParent->obj, popSelectbyVar->obj, popSize1*numObj*sizeof (double));

	// and popSelectbyObj to popParent
	memcpy (popParent->var+popSize1*numVar, popSelectbyObj->var, popSize2*numVar*sizeof (double));
	memcpy (popParent->vel+popSize1*numVar, popSelectbyObj->vel, popSize2*numVar*sizeof (double));
	memcpy (popParent->obj+popSize1*numObj, popSelectbyObj->obj, popSize2*numObj*sizeof (double));

	// CSO
	TSD_reproduce_CSO (pop, popParent, popOffspring); 

	// free popParent
	TSD_Pop_free (popParent);
}

static void TSD_optimize (Population_t* pop) {
	Matrix_t*	Z 	= Population_reference (Problem_get());		// reference vector of objective
	Matrix_t*	V	= NULL;						// reference vector of variable
	int		numObj  = pop->obj->colDim;
	int		numVar  = pop->var->colDim;
	TSD_Pop_t*	popSelectbyVar 	= TSD_Pop_new (Z->rowDim, numVar, numObj);
	TSD_Pop_t*	popSelectbyObj 	= TSD_Pop_new (Z->rowDim, numVar, numObj);
	TSD_Pop_t*	popOffspring 	= TSD_Pop_new (4*Z->rowDim+10, numVar, numObj);
	TSD_EP02_t*	tsdEP02		= NULL;
	TSD_EP03_t*	tsdEP03		= NULL;
	int		EPSize	 	= Z->rowDim;
	int		byVar[numObj+10];
	double		gammaV[Z->rowDim+10]; 
	double		gammaZ[Z->rowDim+10]; 
	int		sel[4*Z->rowDim+10];
	int		i, k;
	
	// if numObj = 2
	if (2 == numObj) { 
		tsdEP02 = TSD_EP02_new (EPSize, numVar, numObj);
	} else { 
		tsdEP03 =  TSD_EP03_new (Z, gammaZ);
	}

	// set byVar 
	for (i=0, k=0; i<numVar; i++) if (1 == pop->CS[i]) byVar[k++] = i;

	// get V
	V = TSD_get_RV (Z->rowDim, k);			

	// set popOffspring size
	popOffspring->popSize = V->rowDim + Z->rowDim;

	// init popOffspring
	TSD_init_pop (pop, popOffspring);

	// get gammaV and gammaZ
	TSD_get_gamma (V, gammaV);
	TSD_get_gamma (Z, gammaZ);

	// set selection flag
	memset (sel, 0, popOffspring->popSize*sizeof (int));

	// selection
	TSD_selection_by_var (popOffspring, V, TSD_get_theta(), gammaV, byVar, sel, popSelectbyVar); 
	TSD_selection_by_obj (popOffspring, Z, TSD_get_theta(), gammaZ, sel, popSelectbyObj);

	// update external population
	if (2 == numObj) {
		TSD_EP02_update (pop, tsdEP02, popSelectbyVar, popSelectbyObj); 
	} else {
		TSD_EP03_update (pop, tsdEP03, popSelectbyVar, popSelectbyObj);
	}

	// loop
	while (0 == isTerminal (pop)) {
		
		// reproduce
		TSD_reproduce (pop, popSelectbyVar, popSelectbyObj, popOffspring);

		// set selection flag
		memset (sel, 0, popOffspring->popSize*sizeof (int));

		// selection
		TSD_selection_by_var (popOffspring, V, TSD_get_theta(), gammaV, byVar, sel, popSelectbyVar);
		TSD_selection_by_obj (popOffspring, Z, TSD_get_theta(), gammaZ, sel, popSelectbyObj);

		// update external population
		if (2 == numObj) {
			TSD_EP02_update (pop, tsdEP02, popSelectbyVar, popSelectbyObj); 
		} else {
			TSD_EP03_update (pop, tsdEP03, popSelectbyVar, popSelectbyObj);
		}
	}
}



/*************************************************************************************************************/
/**************  SLPSO Optimze ********************************************************************************/
/*************************************************************************************************************/

// Aggregation Function
static double g_te (double *f, double *lambda, int len) {
	double 	xValue, t;
	int 	i;
	double*	zIdeal = Problem_getIdealPoint ();
	double	z[len+10];

	// set z
	for (i=0; i<len; i++) if (zIdeal[i] < 0) {
		z[i] = zIdeal[i];
	} else {
		z[i] = 0;
	}

	//
	for (i=0, xValue = 0; i<len; i++) {
		t = lambda[i]*(f[i] - z[i]);
		if (t > xValue)
			xValue = t;
	}

	return xValue;
}

// Terminal condition
static double	SIT_VARIANCE[10] = {1, 2, 3};
static int	SIT_i = 0;

// Terminal condition
static int SLPSO_is_terminal (double *F, int m, int t) {
	int 	i = SIT_i;

	if ((t / m) % 20  == 0) {
		SIT_VARIANCE[i] = VAR (F, m);
		SIT_i = (i+1) % 3;

		if (fabs(2*SIT_VARIANCE[0] - SIT_VARIANCE[1] - SIT_VARIANCE[2]) < DBL_EPSILON) {
			SIT_VARIANCE[0] = 1; SIT_VARIANCE[1] = 2; SIT_VARIANCE[2] = 3;
			return 1;
		}
		if (SIT_VARIANCE[0] < DBL_EPSILON || SIT_VARIANCE[1] < DBL_EPSILON || SIT_VARIANCE[2] < DBL_EPSILON) {
			SIT_VARIANCE[0] = 1; SIT_VARIANCE[1] = 2; SIT_VARIANCE[2] = 3;
			return 1;
		}
	}
	return 0;
}

// SLPSO Optimze
static void SLPSO_optimize (Population_t* pop, double (*g_func)(double *f, double *w, int M), double *gW) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		numVar  = pop->var->colDim;
	int		numObj	= pop->obj->colDim;
	int 		i, j, k, a;
	double		rho = 0.1;
	int		maxFEs = (int)(rho * Problem_getLifetime ());

	//		paramters of SL-PSO
	int 		n = numVar;		// number of variables
	int 		t = 0;
	int 		M = 100; 
	double 		alpha=0.5, belta=0.01;
	int 		m = M + n / 10;		// population size 
	double 		epsilon = belta * n / M;
	double 		PL[m+10]; 
	double 		F[m+10];
	double*		X = (double *)malloc (m*n*sizeof (double));
	double*		V = (double *)calloc (m*n, sizeof (double));
	double*		Y = (double *)calloc (m*numObj, sizeof (double));
	double 		best_solution_x[numVar+10];
	double 		best_solution_y[numObj+10];
	double		Iij, Cij, X_bar[n+10], r1, r2, r3;
	Matrix_t*	T = Matrix_new (m, 1);
	int*		index = NULL;

	// reset population size for 2000, 5000 variables
	m = (m > 200) ? 200 : m; 

	// set best_solution_y
	for (i=0; i<numObj; best_solution_y[i] = 1.0e+100, i++) {};

	// set PL
	for (i=0; i<m; PL[i] = pow(1.0 - 1.0*i/m, alpha*log(ceil(1.0*n/M))), i++) {};

	// init pop 
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++){
			X[i*n+j] = lowBound[j] + randu()*(uppBound[j]-lowBound[j]);
		}
	}

	// main loop
	while (t < maxFEs) {
		// fitness evaluation
		for (i=0; i<m; i++) {
			// make X+i satisfy condition of R
			TSD_satisfy_condition_LS (pop, X+i*numVar);		

			// evaluate 
			Problem_evaluate (X+i*numVar, numVar, Y+i*numObj, numObj);
			t++;
		}

		// aggregation
		for (i=0; i<m; i++) {
			F[i] = g_func (Y+i*numObj, gW, numObj);
		}

		// sort
		memcpy (T->elements, F, m*sizeof (double));
		index = sort (T, (char *)"DES");

		// update best solution using Y
		r1 = g_func (Y+index[m-1]*numObj, gW, numObj);
		r2 = g_func (best_solution_y, gW, numObj);
		if (r1 < r2 ) {
			memcpy (best_solution_x, X+index[m-1]*numVar, numVar*sizeof (double));
			memcpy (best_solution_y, Y+index[m-1]*numObj, numObj*sizeof (double));
			#if (1 == OPT_PRINT)
				printf ("Fitness=%.16e; t=%d/%d; SLPSO\n", r1, t, maxFEs);
			#endif
		}
	
		// is terminal
		if (SLPSO_is_terminal (F, m, t) || t >= maxFEs) { 
			free (index);
			break;	
		}

/***************** update X **************************************************************************************/

		// set X_bar
		mean (X, m, numVar, X_bar);

		// update X_i
		for (a=0; a<m-1; a++) if (randu () < PL[a]) {
			i = index[a];
			for (j=0; j<n; j++) {
				k = rand()%(m-a-1) + (a+1);
				k = index[k];	
				Iij = X[k*n+j] - X[i*n+j];
				Cij = X_bar[j] - X[i*n+j];
				r1 = randu ();
				r2 = randu ();
				r3 = randu ();
				V[i*n+j] = r1*V[i*n+j] + r2*Iij + r3*epsilon*Cij;
				X[i*n+j] = X[i*n+j] + V[i*n+j];
				if (X[i*n+j] > uppBound[j]) {
					X[i*n+j] = uppBound[j];
				}
				if (X[i*n+j] < lowBound[j]) {
					X[i*n+j] = lowBound[j];
				}
			}
		}

		// free index
		free (index);
	}

	// extra the best solution 
	memcpy (pop->var->elements, best_solution_x, numVar*sizeof (double));
	memcpy (pop->obj->elements, best_solution_y, numObj*sizeof (double));

	// free T
	Matrix_free (&T);

	// free X, V, Y
	free (X); free (V); free (Y);
}

/*************************************************************************************************************/
/************* TSD Analysis **********************************************************************************/
/*************************************************************************************************************/

static void TSD_analysis_correlation (Population_t* pop) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		numVar = pop->var->colDim;
	int		numObj = pop->obj->colDim;
	double		X[numVar+10]; 
	double		Y[numObj+10];
	int		i, j;

	// 0. alloc memory for C
	if (!pop->NCOR) {
		pop->NCOR = (int *)calloc(numVar, sizeof (int));
	}

	// 1.  check one by one
	for (i=0; i<numVar; i++) {
		memcpy (X, pop->var->elements, numVar*sizeof (double));
		X[i]  = randu () * (uppBound[i] - lowBound[i]) + lowBound[i];
		Problem_evaluate (X, numVar, Y, numObj);

		// check 
		for (j=0, pop->NCOR[i]=0; j<numObj; j++) {
			if (fabs(Y[j] - pop->obj->elements[j]) > DBL_EPSILON) {
				pop->NCOR[i] += 1;
			}
		}
	}
}

static void TSD_analysis_landscape (Population_t* pop) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		numVar = pop->var->colDim;
	int		numObj = pop->obj->colDim;
	int		i, j, k;

	// 		paramter of perturbation	
	int 		NP = 10;		
	double		U[1000];	// uniform distribution
	double		r1;		// random number
	double		X[numVar+10]; 
	Matrix_t*	Y  = Matrix_new (numVar*NP, numObj); 
	double		L1[(NP+10)*numObj];
	double		L2[(NP+10)*numObj];
	double		L3[(NP+10)*numObj];

	// 0.1 alloc memory for R and C
	if (!pop->R) {
		pop->R = (int *)calloc(numVar*numVar, sizeof (int));
	}

	// 0.2. set one point randomly:	Case one
	for (i=0, r1=randu(); i<numVar; i++) {
		pop->var->elements[i] = lowBound[i] + r1*(uppBound[i] - lowBound[i]);
	}
	Problem_evaluate (pop->var->elements, numVar, pop->obj->elements, numObj);

	// 1 perturb
	for (j=0; j<=NP; U[j] = randu (), j++) {}
	for (i=0; i<numVar; i++) {
		for (j=0; j<NP; j++) {
			memcpy (X, pop->var->elements, numVar*sizeof (double));
			X[i] = lowBound[i] + ((j + U[j])/NP)*(uppBound[i] - lowBound[i]);
			Problem_evaluate (X, numVar, Y->elements+i*NP*numObj+j*numObj, numObj);
		}
	}

	// 2 compute R 
	for (i=0; i<numVar; i++) { 	
		pop->R[i*numVar+i] = 7;
		for (j=i+1; j<numVar; j++) {
			memcpy (L1, Y->elements+i*NP*numObj, NP*numObj*sizeof(double));
			memcpy (L2, Y->elements+j*NP*numObj, NP*numObj*sizeof(double));
			for (k=NP*numObj-1; k>=0; k--) {
				L3[k] = L1[k] - L2[k];			// L3 = L1 - L2
			}
	
			// check if L1 = L2
			if (norm (L3, NP*numObj) < FLT_EPSILON) {
				pop->R[i*numVar+j] = 7;
				pop->R[j*numVar+i] = 7;
			}
		}
	}

	// free Y
	Matrix_free (&Y);
}

static int isConflict (double *A, double *B, int num) {
	int	posi = 0, nega = 0, i;
	double	t;

	for (i=0; i<num; i++) {
		t = A[i] - B[i];
		if (t > DBL_EPSILON) posi++;
		if (t < -DBL_EPSILON) nega++;
	}

	if (posi > 0 && nega > 0) {
		return 1;	
	} 
	return 0;
}

static void TSD_analysis_conflict (Population_t* pop) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		numVar 	= pop->var->colDim;
	int		numObj 	= pop->obj->colDim;
	int		queue[numVar+10], n=0;
	int		i;

	// 		paramter of perturbation	
	double		epsilon = 1.0e-6;		
	double		X[numVar+10]; 
	double 		Y[numObj+10];

	// 0 alloc memory for CV
	if (!pop->CV) {
		pop->CV = (int *)calloc(numVar, sizeof (int));
	} else {
		memset (pop->CV, 0, numVar*sizeof (int));
	} 

	// correlations
	for (i=0, n=0; i<numVar; i++) if (pop->NCOR[i] > 1) {
		queue[n] = i;
		n++;
	}
	if (n < numObj) {
		for (i=0; i<n; i++) {
			pop->CV[queue[i]] = 1;
		}
		return;
	}

	// perturb
	for (i=0; i<numVar; i++) if (pop->NCOR[i] > 1){
		memcpy (X, pop->var->elements, numVar*sizeof (double));
		X[i] = pop->var->elements[i] + epsilon * (uppBound[i] - lowBound[i]);
		if (X[i] < uppBound[i]) {
			Problem_evaluate (X, numVar, Y, numObj);
			if (isConflict(pop->obj->elements, Y, numObj) == 1) {
				pop->CV[i] = 1;
				continue;
			}
		}
		X[i] = pop->var->elements[i] - epsilon * (uppBound[i] - lowBound[i]);
		if (X[i] > lowBound[i]) {
			Problem_evaluate (X, numVar, Y, numObj);
			if (isConflict(pop->obj->elements, Y, numObj) == 1) {
				pop->CV[i] = 1;
				continue;
			}
		}
	}
}

static int CS_isInteractive (Population_t* pop, int i, int j) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		numVar   = pop->var->colDim;
	int		numObj   = pop->obj->colDim;
	double 		x[4*numVar];
	double		y[4*numObj];
	double		a2, b2; 
	double		epsilon = 1.0e-6;
	double 		delta1, delta2;
	int		k;

	memcpy (x, pop->var->elements, numVar*sizeof (double));
	memcpy (y, pop->obj->elements, numObj*sizeof (double));
	a2 = x[i] + epsilon * (uppBound[i] - lowBound[i]);
	if (a2 > uppBound[i]) {
		a2 = x[i] - epsilon * (uppBound[i] - lowBound[i]);
	}
	b2 = x[j] + epsilon * (uppBound[j] - lowBound[j]);
	if (b2 > uppBound[j]) {
		b2 = x[j] - epsilon * (uppBound[j] - lowBound[j]);
	}
	memcpy (x+1*numVar, x, numVar*sizeof (double));
	memcpy (x+2*numVar, x, numVar*sizeof (double));
	memcpy (x+3*numVar, x, numVar*sizeof (double));
	x[1*numVar+i] = a2;
	x[2*numVar+j] = b2;
	x[3*numVar+i] = a2;
	x[3*numVar+j] = b2;
	Problem_evaluate (x+1*numVar, numVar, y+1*numObj, numObj);
	Problem_evaluate (x+2*numVar, numVar, y+2*numObj, numObj);
	Problem_evaluate (x+3*numVar, numVar, y+3*numObj, numObj);
	
	for (k=0; k<numObj; k++) {
		delta1 = y[1*numObj+k] - y[0*numObj+k];
		delta2 = y[3*numObj+k] - y[2*numObj+k];
	
		// base on DG and DG2
		if (fabs(delta1 - delta2) > 1.0e-6) {
			return 1;
		}
	}
	return 0;
}

static void TSD_conflict_space (Population_t* pop) {
	int	numVar  = pop->var->colDim;
	int	numObj  = pop->obj->colDim;
	int 	P[numVar+10], n1=0;	// conflict variables
	int 	Q[numVar+10], n2=0;	// non-conflict variables
	int	I[numVar+10], maxI = 0;
	int 	i, j, a, b, r;
	int	queue[numVar+10], n=0;
	
	//
	if (NULL == pop->CV) return;

	// alloc memory for CS
	if (NULL == pop->CS) {
		pop->CS = (int *)calloc(numVar, sizeof (int));
	} else {
		memset (pop->CS, 0, numVar*sizeof (int));
	}

	// split two cluster: 
	for (i=0; i<numVar; i++) if (pop->CV[i]) {
		P[n1++] = i;
	} else {
		Q[n2++] = i;
	}

	// case one: n1 == 0 or n2 == 0
	if (0 == n1 || 0 == n2) {
		for (i=0, n=0; i<numVar; i++) if (pop->NCOR[i] > 1) {
			queue[n++] = i;
		}
		if (n == 0) {
			r = rand () % numVar;
		} else {
			r = queue[rand () % n];
		}
		pop->CS[r]  = 1;
		return;
	}

	// case two: 0 < n1 < numObj
	if (n1 < numObj) {
		for (i=0; i<n1; i++) { 
			pop->CS[P[i]] = 1;
		}
		return;
	}

	// case three:  n1 >= numObj
	memset (I, 0, numVar*sizeof (int));
	for (a=0; a<n1; a++) {
		i = P[a];
		for (b=0; b<n2; b++) {
			j = Q[b];
			I[i] +=  CS_isInteractive (pop, i, j);
		}
	}
	for (a=0, maxI=0; a<n1; a++) {
		i = P[a];
		if (I[i] > maxI) {
			maxI = I[i];
		}
	}
	for (a=0; a<n1; a++) {
		i = P[a];
		if (I[i] == maxI) {
			pop->CS[i] = 1;
		}
	}
}

/*************************************************************************************************************/
/*************	Satisfy Conditions ***************************************************************************/
/*************************************************************************************************************/
// 
static int SCL_index[1000000] = {-1};
//
static void TSD_satisfy_condition_LS (Population_t*pop, double *x) {
	double*	lowBound = Problem_getLowerBound ();
	double*	uppBound = Problem_getUpperBound ();
	int	numVar  = pop->var->colDim;
	int 	i, j;

	if (NULL == pop->R)  return;
	if (NULL == pop->CV) return;

	if (-1 == SCL_index[0]) for (i=numVar-1; i>=0; i--) {
		for (j=0; j<=i; j++) if (7 == pop->R[i*numVar+j]) {
			SCL_index[i] = j;
			break;
		}
	}

	for (i=0; i<numVar; i++) if (0 == pop->CV[i] && i != SCL_index[i]) {
		j = SCL_index[i];
		x[i] = ((x[j] - lowBound[j])/(uppBound[j] - lowBound[j]))*(uppBound[i] - lowBound[i]) + lowBound[i];
	}
}

static void TSD_satisfy_condition_CS (Population_t*pop, double *parent, double* child) {
	int	numVar  = pop->var->colDim;
	int	i, r = rand () % 100;
	double	d;
	
	if (NULL == pop->CS) return;

	if (r < 5) {
		for (i=0, d=0; i<numVar; i++) if (1 == pop->CS[i]) {
			d =+ fabs(parent[i] - child[i]);
			if (d > DBL_EPSILON) break;
		}
		if (d > DBL_EPSILON) {
			for (i=0; i<numVar; i++) if (0 == pop->CS[i]) {
				child[i] = parent[i];
			}
		}
	} else if (r < 50) {
		for (i=0, d=0; i<numVar; i++) if (0 == pop->CS[i]) {
			d =+ fabs(parent[i] - child[i]);
			if (d > DBL_EPSILON) break;
		}
		if (d > DBL_EPSILON) {
			for (i=0; i<numVar; i++) if (1 == pop->CS[i]) {
				child[i] = parent[i];
			}
		}
	}
}

/*************************************************************************************************************/
/************************** END ******************************************************************************/
/*************************************************************************************************************/
