#ifndef _BPMLP_H
#define _BPMLP_H

#include "matrix.h"

typedef struct BPMLP_tag {
	double**	layer;		// layer
	double** 	error;	
	int		layer_dim1;
	int*		layer_dim2;

	double***	weight;		// weight
	int		weight_dim1;
	int*		weight_dim2;
	int** 		weight_dim3;

	double 	rate;	// learning rate	

} BPMLP_t;

BPMLP_t*	BPMLP_new (int inputNeuron, int middleNeuron, int outputNeuron, int layerNum, double rate);
void 		BPMLP_free (BPMLP_t** pbpmlp);
double*		BPMLP_computeOutput (BPMLP_t* bpmlp, double* input);
void 		BPMLP_training (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, double p, int epochs);
void 		BPMLP_training (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, int epochs);
double		BPMLP_training_by_PSO (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target, int maxFEs);
double		BPMLP_evaluate (BPMLP_t* bpmlp, Matrix_t* input, Matrix_t *target);

#endif 
