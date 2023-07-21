#include "model.h"
#include "population.h"
#include "problem.h"
#include "myrandom.h"
#include "algebra.h"
#include "bpmlp.h"
#include <math.h>
#include <string.h>
#include <float.h>

#define MODEL_EPSILON 1.0e-4
#define MODEL_PRINT 0

/*************************************************************************************************************/
/**************  Basic function  *****************************************************************************/
/*************************************************************************************************************/

// linear function
static double func100 (double *x, int n, double *para) {
	return para[n];
}

static double func101 (double *x, int n, double *para) {
	double 	y = 0;	
	int 	i;

	for (i=0; i<n; i++) { 
		y += x[i]*para[i]; 
	}
	return y;
}

static double func102 (double *x, int n, double *para) {
	double 	y = 0;	
	int 	i;

	for (i=0; i<n; i++) { 
		y += x[i]*para[i]; 
	}
	y += para[n];
	return y;
}

/*************************************************************************************************************/

// segment function ||
static double func200 (double *x, int n, double *para) { 
	double 	y;
	
	y = func102 (x, n, para); 
	y = fabs (y);
	return  y;
}

/*************************************************************************************************************/

// sin (x)
static double func300 (double *x, int n, double *para) { 
	double 	y;
	
	y = func102 (x, n, para); 
	y = 0.5*sin (y) + 0.5;
	return y; 
}

/*************************************************************************************************************/

// y = x^n
static double func400 (double *x, int n, double *para) { 
	double 	y;
	
	y = func102 (x, n, para); 
	y = fabs (y);
	y = pow (y, para[n+1]); 
	return y;
}

/*************************************************************************************************************/

static double (* funcArr[1000])(double *x, int n, double *para) = {func100, func101, func102, func200, func300, func400};
int 	numFun = 3+1+1+1;	

/*************************************************************************************************************/
/************** evaluate parameter of Model ******************************************************************/
/*************************************************************************************************************/

static double evaluate (double *para, Matrix_t *X, double *Y, double (*func)(double*,int,double*)) {
	int 	rowDim = X->rowDim;
	int 	colDim = X->colDim;
	double 	MSE = 0, y;
	int 	i;

	for (i=0; i<rowDim; i++) {
		y = func (X->elements+i*colDim, colDim, para);
		MSE += (y - Y[i])*(y - Y[i]);
	}
	if (rowDim > 0) {
		MSE /= rowDim;
	}
	return MSE;
}

/*************************************************************************************************************/
/************** Model ****************************************************************************************/
/*************************************************************************************************************/

typedef struct Model_tag {
	int		type; 		// type of model, 0: basic function; 2: MLP
	double		(*func) (double *x, int n, double* para);
	double*		para;	
	BPMLP_t*	mlp;		// MLP reuron Netwrok
	double 		MSE;		// Mean-squareerror (MSE)
} Model_t;

static Model_t* model_pool;
static double	para_lowBound[1000];
static double	para_uppBound[1000];

/*************************************************************************************************************/
/************** utilization tools ****************************************************************************/
/*************************************************************************************************************/

static char *toString (double (*func)(double*, int, double*)) {
	if (func100 == func) {
		return (char *)"func100";
	} else if (func101 == func) {
		return (char *)"func101";
	} else if (func102 == func) {
		return (char *)"func102";
	} else if (func200 == func) {
		return (char *)"func200";
	} else if (func300 == func) {
		return (char *)"func300";
	} else if (func400 == func) {
		return (char *)"func400";
	} else {
		return (char *)"func-xxx";
	}
}

void model_print (Model_t* model) {
	int 	i;
	
	if (0 == model->type) {
		printf ("Model.type = Basic Function\n");
		printf ("Model.func = %s\n", toString (model->func));
		for (i=0; i<10; i++) {
			printf ("Para[%d] = %.16e\n", i, model->para[i]);
		}
	} else {
		printf ("Model.type = MLP\n");
	}
	printf ("Model.MSE = %.16e\n", model->MSE);
}

static int IR_flag[10000000];
static int IR_len; 
static int is_remove (Population_t* pop, int cursor, int depV) {
	double*	lowBound = Problem_getLowerBound ();
	double*	uppBound = Problem_getUpperBound ();
	int	numVar = pop->var->colDim;
	int	numObj = pop->obj->colDim;
	double 	x[numVar+10];
	double 	y[numObj+10];
	int 	i, j;

	//
	for (i=IR_len; i<=cursor; i++) {
		for (j=0; j<numVar; j++) {
			memcpy (x, pop->var->elements+i*numVar, numVar*sizeof (double));
			x[j] = lowBound[j] + randu()*(uppBound[j] - lowBound[j]);
			Problem_evaluate (x, numVar, y, numObj);
			if (distance_p2p (y, pop->obj->elements+i*numObj, numObj) < FLT_EPSILON) {
				IR_flag[i*numVar+j] = 1;
			}
		}
	}
	if (cursor+1 > IR_len) {
		IR_len = cursor + 1;
	}

	return IR_flag[cursor*numVar+depV];	
}

/*************************************************************************************************************/
/********** Set Paramter Boundary  ***************************************************************************/
/*************************************************************************************************************/

static void set_para_boundary (double (*func)(double*, int, double*), int numVar) {
	int 	i;

	for (i=0; i<numVar+2; i++) {
		para_uppBound[i] =  100;
		para_lowBound[i] = -100;
	}
	
	if (func == func200) {			// |x|
		para_uppBound[numVar] =  2;
		para_lowBound[numVar] = -2;
	} else if (func == func300) {		// sin (x)
		para_uppBound[numVar] = 3.15; 
		para_lowBound[numVar] = 0;
	} else if (func == func400) { 		// x^n
		para_uppBound[numVar] =  2;
		para_lowBound[numVar] = -2;
		para_uppBound[numVar+1] = 100;
		para_lowBound[numVar+1] = 0;
	}
#if (1 == MODEL_PRINT)
	for (i=0; i<numVar+2; i++) {
		printf("%s: (low, upp)[%d] = %f, %f\n", toString(func), i, para_lowBound[i], para_uppBound[i]);
	}
#endif
}

/*************************************************************************************************************/
/*************** get parameters  *****************************************************************************/
/*************************************************************************************************************/

static double	SIT_VARIANCE[10] = {1, 2, 3};
static int	SIT_i = 0;
static int SLPSO_is_terminal (double* F, int m, int t) {
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

static void when_func100 (double *para, double &MSE, Matrix_t* DATA, double *R) {
	int 	colDim = DATA->colDim;
	int 	i;
	
	for (i=0; i<colDim; i++) {
		para[i] = 0;
	}
	para[colDim] = R[0];
	MSE = evaluate (para, DATA, R, func100); 
}

static void when_func101 (double *para, double &MSE, Matrix_t* DATA, double *R) {
	int 	colDim = DATA->colDim;
	int 	rowDim = DATA->rowDim;
	double	x[rowDim+10], A, B;
	int 	i, j;
	
	B = mean (R, rowDim);
	for (i=0; i<colDim; i++) {
		para[i] = 0;
	}
	for (i=0; i<colDim; i++) {
		for (j=0; j<rowDim; j++) {
			x[j] = DATA->elements[j*colDim+i];
		}
		A = mean (x, rowDim);
		if (fabs(A) > DBL_EPSILON) {
			para[i] = B / A;
			break;
		}
	}
	para[colDim] = 0;
	MSE = evaluate (para, DATA, R, func101); 
}

static void when_func102 (double *para, double &MSE, Matrix_t* DATA, double *R) {
	int		rowDim = DATA->rowDim;
	int 		colDim = DATA->colDim;
	Matrix_t*	XI = Matrix_new (rowDim, colDim+1);
	double  	belta[colDim+10], delta;
	int 		i, j;

	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			XI->elements[i*(colDim+1)+j] = DATA->elements[i*colDim+j];
		}
		XI->elements[i*(colDim+1) + j] = 1.0;
	}

	// make line
	make_line (XI, R, belta, delta);

	// set para
	for (i=0; i<=colDim; i++) {
		para[i] = belta[i];
	}
	MSE = evaluate (para, DATA, R, func102); 

	// free XI
	Matrix_free (&XI);
}

static void 
get_func_para_by_SLPSO (double *para, double &MSE, Matrix_t* DATA, double *R, double (*func)(double*, int, double*)) {
	int 		i, j, k, a, b;

	//	paramters of SL-PSO
	int 		n = DATA->colDim + 2;	// number of variable
	int 		t = 0;
	int 		M = 100; 
	double 		alpha=0.5, belta=0.01;
	int 		m = M + n / 10; 	// population size
	double 		epsilon = belta * n / M;
	double 		PL[m+10]; 
	double		F[m+10];
	double 		best_solution_x[n+10];
	double 		best_solution_y;
	double*		X = (double *)malloc (m*n*sizeof (double));
	double*		V = (double *)calloc (m*n, sizeof (double));
	double		Iij, Cij, X_bar[n+10], r1, r2, r3;
	Matrix_t*	T = Matrix_new (m, 1);
	int*		index = NULL;
	int 		maxFEs = 1000000;

	// deal with linear function
	if (func100 == func) {
		when_func100 (para, MSE, DATA, R);
	} else if (func101 == func) {
		when_func101 (para, MSE, DATA, R);
	} else if (func102 == func) {
		when_func102 (para, MSE, DATA, R);
	}
	if (func100 == func || func101 == func || func102 == func) {
		Matrix_free (&T); free (X); free (V);
		return;
	}

	// set PL
	for (i=0; i<m; i++) { 
		PL[i] = pow(1.0 - 1.0*i/m, alpha*log(ceil(1.0*n/M))); 
	}

	// set best solution
	best_solution_y = 1.0e+100;

	// set parameter boundary 	
	set_para_boundary (func, DATA->colDim);

	// init pop
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			X[i*n+j] = para_lowBound[j] + randu()*(para_uppBound[j] - para_lowBound[j]);
		}
	}

	// main loop
	while (t < maxFEs) {

		// fitness evaluation
		for (i=0; i<m; i++) { 
			F[i] = evaluate (X+i*n, DATA, R, func); 
			t++;
		}

		// set X_bar
		mean (X, m, n, X_bar);

		// sort
		memcpy (T->elements, F, m*sizeof (double));
		index = sort (T, (char *)"DES");

		// update best solution 
		if (F[index[m-1]] < best_solution_y) {
			memcpy (best_solution_x, X+index[m-1]*n, n*sizeof (double));
			best_solution_y = F[index[m-1]];

		#if (1 == MODEL_PRINT)
			printf ("%s; F = %e; %d/%d\n", toString(func), best_solution_y, t, maxFEs);
		#endif
		}
		
		if (SLPSO_is_terminal (F, m, t)) {
			free (index);
			break;
		}
		
		// update X_i
		for (a=0; a<m-1; a++) if (randu () < PL[a]) {
			i = index[a];
			for (j=0; j<n; j++) {
				b = rand()%(m-a-1) + (a+1);
				k = index[b];	
				Iij = X[k*n+j] - X[i*n+j];
				Cij = X_bar[j] - X[i*n+j];
				r1 = randu ();
				r2 = randu ();
				r3 = randu ();
				V[i*n+j] = r1*V[i*n+j] + r2*Iij + r3*epsilon*Cij;
				X[i*n+j] = X[i*n+j] + V[i*n+j];
				if (X[i*n+j] < para_lowBound[j]) {
					 X[i*n+j] = para_lowBound[j]; 
				}
				if (X[i*n+j] > para_uppBound[j]) { 
					X[i*n+j] = para_uppBound[j]; 
				}
			}

		}
		// free index
		free (index);
	}

	// extra best solution
	memcpy (para, best_solution_x, n*sizeof (double));
	MSE = best_solution_y;

	// free T
	Matrix_free (&T); 

	// free X, V
	free (X); free (V);
}

/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

// model update
void model_update (Model_t* model, Matrix_t* X, double* Y) {
	int		rowDim = X->rowDim;	
	int		colDim = X->colDim;	
	double 		para[colDim+30];
	double 		MSE = 0;
	int 		i, j, k;
	double		t;
	BPMLP_t*	mlp = NULL;
	Matrix_t*	mlp_input = NULL;	
	Matrix_t*	mlp_targt = NULL;	

	// allocate memory
	if (NULL == model->para) {
		model->para = (double *)calloc (colDim+20, sizeof (double));
	}

	// if no data
	if (NULL == X || X->rowDim < 3) return;

	// set para
	model->MSE = 1.0e+100;
	for (i=0; i<numFun; i++) {
		for (j=0; j<=rowDim; j++) if (model->MSE > FLT_EPSILON) {

			// remove j-th row for i = 0, 1, 2
			if (j<rowDim && i<3) {
				for (k=0; k<colDim; k++) {
					t = X->elements[j*colDim+k];
					X->elements[j*colDim+k] =  X->elements[(rowDim-1)*colDim+k];
					X->elements[(rowDim-1)*colDim+k] = t;
				}
				t = Y[j];
				Y[j] = Y[rowDim-1];
				Y[rowDim-1] = t;
				X->rowDim -= 1; 	// X->rowDim--;
			}

			// training
			memset (para, 0, (colDim+10)*sizeof (double));
			get_func_para_by_SLPSO (para, MSE, X, Y, funcArr[i]);
			if (MSE < model->MSE) {
				model->type = 0;	
				model->func = funcArr[i];
				memcpy (model->para, para, (colDim+10)*sizeof (double));
				model->MSE = MSE;
			}

			// recover j-th row for i = 0, 1, 2
			if (j<rowDim && i<3) {
				for (k=0; k<colDim; k++) {
					t = X->elements[j*colDim+k];
					X->elements[j*colDim+k] =  X->elements[(rowDim-1)*colDim+k];
					X->elements[(rowDim-1)*colDim+k] = t;
				}
				t = Y[j];
				Y[j] = Y[rowDim-1];
				Y[rowDim-1] = t;
				X->rowDim += 1; 	// X->rowDim++;
			}
		}
	}

	if (model->MSE > MODEL_EPSILON) {
		if (NULL == model->mlp) {
			model->mlp = BPMLP_new (X->colDim, 10, 1, 3, 0.01);	 
		}

		// training 	
		mlp = model->mlp;
		mlp_input = X;
		mlp_targt = Matrix_new (X->rowDim, 1);	
		memcpy (mlp_targt->elements, Y, X->rowDim*sizeof (double));
		MSE = BPMLP_training_by_PSO (mlp, mlp_input, mlp_targt, 1000000);
		if (MSE < model->MSE) {
			model->type = 1;	
			model->MSE = MSE;
		}
		Matrix_free (&mlp_targt);
	}
}

/*************************************************************************************************************/
/*************************************************************************************************************/
/*************************************************************************************************************/

// model update
void model_update (Population_t *pop) {
        int    		numVar  = pop->var->colDim;
        int    		popSize = pop->cursor;
	Matrix_t*	X = NULL; 
	double		Y[popSize+10];
	int		i, j, k, n, v;
	int 		argc = 0;

	// normalize decision variables
	Problem_xToOne (pop->var);
	
	// allocate memory 
	if (NULL == model_pool) {
		model_pool = (Model_t *)calloc (numVar+10, sizeof (Model_t));
	}
	
	// compute argc
	for (v=0, argc=0; v<numVar; v++) if (0 == pop->I[v]) argc++;

	//
	X = Matrix_new (popSize+10, argc);	

	// 1. update model
	for (v=0; v<numVar; v++) {
		// get X, Y
		for (i=0, n=0; i<popSize; i++) if (0 == is_remove(pop, i, v)) {
			// set X
			for (j=0, k=0; j<numVar; j++) if (0 == pop->I[j]) {
				X->elements[n*argc+k] = pop->var->elements[i*numVar+j];
				k++;
			}

			// set Y
			Y[n] = pop->var->elements[i*numVar+v];
			n++;
			if (n >= 20) break;
		}
		X->rowDim = n;

	#if (1 == MODEL_PRINT)
		printf ("%d: X | Y\n", v);
		for (i=0; i<X->rowDim; i++) {
			for (j=0; j<X->colDim; j++) {
				printf ("%f ", X->elements[i*X->colDim+j]);
			}
			printf ("| %f\n", Y[i]);
		}
	#endif

		// update model
		model_update (&model_pool[v], X, Y);

	#if (1 == MODEL_PRINT)
		model_print (&model_pool[v]);
	#endif
	}

	// back from normalization
	Problem_xFromOne (pop->var);
}

/*************************************************************************************************************/
/********** Model get Value **********************************************************************************/
/*************************************************************************************************************/

double model_getValue (Population_t *pop, int depV, double* offset) {
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	int		numVar = pop->var->colDim;
	double 		x[numVar];
	int 		argc = 0;
	int 		j, k;
	double*		mlp_output = NULL;
	double 		value;


	if (NULL == model_pool ||  NULL == model_pool[depV].func) {
		value  = (uppBound[depV] + lowBound[depV]) / 2; 
		if (offset != NULL) {
			offset[0] = (uppBound[depV] - lowBound[depV]) / 2; 
		}
		return value;
	}

	// normalize
	Problem_xToOne (pop->var->elements+pop->cursor*numVar, numVar);

	// set x
	for (j=0, k=0; j<numVar; j++) if (0 == pop->I[j]) {
		x[k] = pop->var->elements[pop->cursor*numVar+j];
		k++;
	}
	argc = k;

	// get value
	if (0 == model_pool[depV].type) {
		value = model_pool[depV].func(x, argc, model_pool[depV].para);
	} else {
		mlp_output = BPMLP_computeOutput (model_pool[depV].mlp, x);
		value = mlp_output[0];
	}
	if (offset != NULL) {
		offset[0]= sqrt (model_pool[depV].MSE);
	}

	// back from normalization
	Problem_xFromOne (pop->var->elements+pop->cursor*numVar, numVar);
	value = value * (uppBound[depV] - lowBound[depV]) + lowBound[depV];
	return value;
}
