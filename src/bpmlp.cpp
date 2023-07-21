#include "bpmlp.h"
#include "myrandom.h"
#include "matrix.h"
#include "algebra.h"
#include <float.h>
#include <string.h>
#include <math.h>

// #define PRINT_WEIGHT 

BPMLP_t* BPMLP_new (int inputNeuron, int middleNeuron, int outputNeuron, int layerNum, double rate) {
	BPMLP_t*	bpmlp = (BPMLP_t *)calloc (1, sizeof (BPMLP_t));
	int		i, j, k;

	if (layerNum % 2 == 0) layerNum++;	// The total number of layers of the auto-encoder is odd

	bpmlp->layer_dim1 = layerNum;
	bpmlp->layer_dim2 = (int *)calloc (bpmlp->layer_dim1, sizeof (int));
	bpmlp->layer_dim2[0] = inputNeuron;
	bpmlp->layer_dim2[layerNum/2] = middleNeuron;
	bpmlp->layer_dim2[layerNum-1] = outputNeuron;

	for (i=1; i<layerNum/2; i++) {
		bpmlp->layer_dim2[i] = bpmlp->layer_dim2[i-1] - (inputNeuron - middleNeuron) / (layerNum/2);
	}
	for (i=layerNum-2; i>layerNum/2; i--) {
		bpmlp->layer_dim2[i] = bpmlp->layer_dim2[i+1] - (outputNeuron - middleNeuron) / (layerNum/2);
	}

#ifdef PRINT_WEIGHT
	// for test
	for (i=0; i<bpmlp->layer_dim1; i++) {
		printf ("%d : %d\n", i, bpmlp->layer_dim2[i]);
	}
#endif

	// allocate memory for layer
	bpmlp->layer = (double **) calloc (bpmlp->layer_dim1, sizeof (double *));
	bpmlp->error = (double **) calloc (bpmlp->layer_dim1, sizeof (double *));
	for (i=0; i<bpmlp->layer_dim1; i++) {
		bpmlp->layer[i] = (double *) calloc (bpmlp->layer_dim2[i], sizeof (double));
		bpmlp->error[i] = (double *) calloc (bpmlp->layer_dim2[i], sizeof (double));
	}
	

	bpmlp->weight_dim1 = layerNum-1;
	bpmlp->weight_dim2 = (int *)calloc (bpmlp->weight_dim1, sizeof (int));
	for (i=0; i<bpmlp->weight_dim1; i++) {
		bpmlp->weight_dim2[i] = bpmlp->layer_dim2[i] + 1;
	}
	bpmlp->weight_dim3 = (int **)calloc (bpmlp->weight_dim1, sizeof (int *));
	for (i=0; i<bpmlp->weight_dim1; i++) {
		bpmlp->weight_dim3[i] = (int *)calloc (bpmlp->weight_dim2[i], sizeof (int));
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			bpmlp->weight_dim3[i][j] = bpmlp->layer_dim2[i+1] + 1;
		}
	}

#ifdef PRINT_WEIGHT
	// for test
	for (i=0; i<bpmlp->weight_dim1; i++) {
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			printf ("%d %d : %d\n", i, j, bpmlp->weight_dim3[i][j]);
		}
	}
#endif

	// allocate memory for weight
	bpmlp->weight = (double ***)calloc (bpmlp->weight_dim1, sizeof (double **));
	for (i=0; i<bpmlp->weight_dim1; i++) {
		bpmlp->weight[i] = (double **)calloc (bpmlp->weight_dim2[i], sizeof (double *));
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			bpmlp->weight[i][j] = (double *)calloc (bpmlp->weight_dim3[i][j], sizeof (double));
			for (k=0; k<bpmlp->weight_dim3[i][j]; k++) {
				bpmlp->weight[i][j][k] = randu ();
#ifdef PRINT_WEIGHT
				printf ("%d, %d, %d: %f\n", i, j, k, bpmlp->weight[i][j][k]);
#endif
			}
		}
	}

	// set learning rate
	bpmlp->rate = rate;

	return bpmlp;	
}

void BPMLP_free (BPMLP_t** pbpmlp) {
	int	 i, j;
	BPMLP_t* bpmlp = *pbpmlp;

	for (i=0; i<bpmlp->weight_dim1; i++) {
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			free (bpmlp->weight[i][j]);
		}
		free (bpmlp->weight[i]);
	}
	free (bpmlp->weight);

	for (i=0; i<bpmlp->weight_dim1; i++) {
		free (bpmlp->weight_dim3[i]);
	}
	free (bpmlp->weight_dim3);
	free (bpmlp->weight_dim2);

	for (i=0; i<bpmlp->layer_dim1; i++) {
		free (bpmlp->layer[i]);
		free (bpmlp->error[i]);
	}
	free (bpmlp->layer);
	free (bpmlp->error);
	free (bpmlp->layer_dim2);
	free (bpmlp);
}

static void resetWeights (BPMLP_t* bpmlp) {
	int i, j, k;

	for (i=0; i<bpmlp->weight_dim1; i++) {
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			for (k=0; k<bpmlp->weight_dim3[i][j]; k++) {
				bpmlp->weight[i][j][k] = randu ();
			}
		}
	}
}

static double sigmoid(double z) {
	 return 1.0 / ( 1 + exp(-z));
 }

double* BPMLP_computeOutput (BPMLP_t* bpmlp, double* input) {
	int 		i, j, k;
	double 		z;
	double*** 	weight = bpmlp->weight;
	double**  	layer  = bpmlp->layer; 

	for (i=0; i<bpmlp->layer_dim2[0]; i++) {
		layer[0][i] = input[i];
	}

	for(i=1; i<bpmlp->layer_dim1; i++) {				//start from the first hidden layer, i.e., i = 1
		for(k=0; k<bpmlp->layer_dim2[i]; k++) {
                   	z = weight[i-1][bpmlp->layer_dim2[i-1]][k];	//z is initialized as the bias
                       	for(j=0; j<bpmlp->layer_dim2[i-1]; j++) {
                              	z += weight[i-1][j][k]*layer[i-1][j]; 	//get the weighted sum
                       	}
                       	layer[i][k] = sigmoid(z);			// active function
              	}
      	}

        return (double *)layer[bpmlp->layer_dim1-1];
}

static void backPropagation(BPMLP_t* bpmlp, double* target) {
	int 		i, j, k;
        double 		err = 0.0;
	double**	layer = bpmlp->layer;
	double**	error = bpmlp->error;
	double***	weight = bpmlp->weight;
	double		rate = bpmlp->rate;

	//Calculate the error of the last layer first
      	i = bpmlp->layer_dim1 - 1;	//the last layer, i.e., the output layer
        for(j=0; j<bpmlp->layer_dim2[i]; j++) {
		//here use (target[j]-layer[i][j]), so the new weights is updates by adding rate*error*layer 
		error[i][j] = layer[i][j]*(1-layer[i][j])*(target[j]-layer[i][j]); 
	}

        while(i-- > 0){	//until to the fist layer, i.e., the input layer
		//The error and the weight are updated simultaneously
                for(j=0; j<bpmlp->layer_dim2[i]; j++) {
                       	err = 0.0;
                      	for(k=0; k<bpmlp->layer_dim2[i+1]; k++) {
                             	err += weight[i][j][k]*error[i+1][k];
                                weight[i][j][k] += rate*error[i+1][k]*layer[i][j];
                               	if(j == bpmlp->layer_dim2[i]-1) {	//Adjust the bias
                                	weight[i][j+1][k] += rate*error[i+1][k];
				}
                        }
                        error[i][j] = err*layer[i][j]*(1.0-layer[i][j]);
               	}
	}
}

//If the mean squared error (MSE) is less than p, it's correct
static int compare (double* a1, double* a2, double p, int len) {
	int 	i;
 	double 	err = 0.0;

 	for(i=0; i<len; i++) {
	 	err += (a1[i]-a2[i])*(a1[i]-a2[i]);
	}
 	err /= len;
 	if(err < p) return 1;
	return 0;
}
 
static double reportModel(BPMLP_t* bpmlp, Matrix_t *input , double p) {
	int 	i, count=0;
	double 	*result = NULL;
	
	for(i=0; i<input->rowDim; i++) {
		result = BPMLP_computeOutput(bpmlp, input->elements+i*input->colDim);
 		if( compare(result, input->elements+i*input->colDim, p, input->colDim)) ++count;
	}
	return count/(double)input->rowDim;
}

void BPMLP_training (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, double p, int epochs) {
	int 	xn=0, i;
        double 	tmp;

       	do {
		for(i=0; i<input->rowDim; i++) {
                	//Forward passing process
                     	BPMLP_computeOutput(bpmlp, input->elements+i*input->colDim);

                     	//Backward adjusting process
                     	backPropagation(bpmlp, target->elements+i*target->colDim);
                }
              	tmp = reportModel(bpmlp, input, p);
                ++xn;
                if(xn % 100 == 0) {
	                printf("Iterative training %d times. The current accuracy of the model is：%f\n", xn, tmp);
                    	p += 0.01; 	//Progressively increase the value of p,
               	}
	} while (tmp < 1.0-p && xn < epochs);
}

void BPMLP_training (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, int epochs) {
	int t, i;

 	resetWeights(bpmlp);
	for(t=0; t<epochs; t++) {
		 for(i=0; i<input->rowDim; i++) {
			 //Forward passing process
			 BPMLP_computeOutput(bpmlp, input->elements+i*input->colDim);

			 //Backward adjusting process
			 backPropagation(bpmlp, target->elements+i*target->colDim);
		}
	}
}

static double evaluate (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, double *x) {
	int 	i, j, k;
	int	num = 0;
	double*	output = NULL;
	double 	MSE = 0, d=0;
	
	for (i=0, num=0; i<bpmlp->weight_dim1; i++) {
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			for (k=0; k<bpmlp->weight_dim3[i][j]; k++) {
				bpmlp->weight[i][j][k] = x[num++];
			}
		}
	}
	for (i=0; i<input->rowDim; i++) {
		output = BPMLP_computeOutput (bpmlp, input->elements+i*input->colDim);
		for (j=0; j<target->colDim; j++) {
			d = output[j] - target->elements[i*target->colDim+j];
			MSE += d * d;
		}
	}
	MSE /= input->rowDim * target->colDim;
	
	return MSE;
}

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

double BPMLP_training_by_PSO (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, int maxFEs) {
	int 	i, j, k, a, b;
	int	num = 0;	

	// compute number of weight
	for (i=0, num=0; i<bpmlp->weight_dim1; i++) {
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			num += bpmlp->weight_dim3[i][j];
		}
	}

	//		paramters of SL-PSO
	int 		n = num;		// number of variable
	int 		t = 0;
	int 		M = 100; 
	double 		alpha=0.5, belta=0.01;
	int 		m = M; 			// population size
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
	double		lowBound[n+10];
	double		uppBound[n+10];

	// set PL
	for (i=0; i<m; i++) { 
		PL[i] = pow(1.0 - 1.0*i/m, alpha*log(ceil(1.0*n/M))); 
	}

	// set best solution
	best_solution_y = 1.0e+100;

	// set boundary
	for (i=0; i<n; i++) {
		uppBound[i] =  10;
		lowBound[i] = -10;
	}

	// init pop
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			X[i*n+j] = lowBound[j] + randu()*(uppBound[j] - lowBound[j]);
		}
	}

	// main loop
	while (t < maxFEs) {

		// fitness evaluation
		for (i=0; i<m; i++) { 
			F[i] = evaluate (bpmlp, input, target, X+i*n);
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
		}
		
		if (SLPSO_is_terminal (F, m, t) || t >= maxFEs) {
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
				if (X[i*n+j] < lowBound[j]) {
					 X[i*n+j] = lowBound[j]; 
				}
				if (X[i*n+j] > uppBound[j]) { 
					X[i*n+j] = uppBound[j]; 
				}
			}

		}
		// free index
		free (index);
	}

	// extra best solution
	for (i=0, num=0; i<bpmlp->weight_dim1; i++) {
		for (j=0; j<bpmlp->weight_dim2[i]; j++) {
			for (k=0; k<bpmlp->weight_dim3[i][j]; k++) {
				bpmlp->weight[i][j][k] = best_solution_x[num++];
			}
		}
	}

	// free T
	Matrix_free (&T); 

	// free X, V
	free (X); free (V);

	return best_solution_y;
}


double BPMLP_evaluate (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target) {
	int 	i, j;
	double*	output = NULL;
	double 	MSE = 0, d=0;
	
	for (i=0; i<input->rowDim; i++) {
		output = BPMLP_computeOutput (bpmlp, input->elements+i*input->colDim);
		for (j=0; j<target->colDim; j++) {
			d = output[j] - target->elements[i*target->colDim+j];
			MSE += d * d;
		}
	}
	if (input->rowDim * target->colDim > 0) {
		MSE /= input->rowDim * target->colDim;
	}
	
	return MSE;
}
