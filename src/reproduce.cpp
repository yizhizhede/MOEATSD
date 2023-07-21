#include "reproduce.h"
#include "matrix.h"
#include "recombination.h"
#include "myrandom.h"
#include "tournament.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

Population_t *Reproduce_one_by_random (Population_t *P) {
	Population_t *Q = (Population_t *)malloc (sizeof (Population_t));
	Q->var = Matrix_new (P->var->rowDim+1, P->var->colDim);
	Q->obj = Matrix_new (P->obj->rowDim+1, P->obj->colDim);

	SBX_reproduce (P->var->elements + ((int)(P->var->rowDim*randu()))*P->var->colDim, 
		       P->var->elements + ((int)(P->var->rowDim*randu()))*P->var->colDim, 
		       Q->var->elements, 
		       Q->obj->elements);

	memcpy (Q->var->elements+Q->var->colDim, P->var->elements, P->var->rowDim*P->var->colDim*sizeof (double));
	memcpy (Q->obj->elements+Q->obj->colDim, P->obj->elements, P->obj->rowDim*P->obj->colDim*sizeof (double));
		
	return Q;
}

Population_t *Reproduce_equ_by_random (Population_t *P) {
	int i;
	Population_t *Q = (Population_t *)malloc (sizeof (Population_t));
	Q->var = Matrix_new (2*P->var->rowDim, P->var->colDim);
	Q->obj = Matrix_new (2*P->obj->rowDim, P->obj->colDim);

	for (i=0; i<P->obj->rowDim; i++) {
		SBX_reproduce (P->var->elements + ((int)(P->var->rowDim*randu()))*P->var->colDim, 
			       P->var->elements + ((int)(P->var->rowDim*randu()))*P->var->colDim, 
			       Q->var->elements + i * Q->var->colDim, 
			       Q->obj->elements + i * Q->obj->colDim);
	}

	memcpy (Q->var->elements + i*Q->var->colDim, P->var->elements, P->var->rowDim*P->var->colDim*sizeof (double));
	memcpy (Q->obj->elements + i*Q->obj->colDim, P->obj->elements, P->obj->rowDim*P->obj->colDim*sizeof (double));
		
	return Q;
}

Population_t *Reproduce_equ_by_tournament (Population_t *P) {
	int i;
	Population_t *Q = (Population_t *)malloc (sizeof (Population_t));
	Q->var = Matrix_new (2*P->var->rowDim, P->var->colDim);
	Q->obj = Matrix_new (2*P->obj->rowDim, P->obj->colDim);

	for (i=0; i<P->obj->rowDim; i++) {
		SBX_reproduce (P->var->elements + tournament(P->var->rowDim)*P->var->colDim, 
			       P->var->elements + tournament(P->var->rowDim)*P->var->colDim, 
			       Q->var->elements + i * Q->var->colDim, 
			       Q->obj->elements + i * Q->obj->colDim);
	}

	memcpy (Q->var->elements + i*Q->var->colDim, P->var->elements, P->var->rowDim*P->var->colDim*sizeof (double));
	memcpy (Q->obj->elements + i*Q->obj->colDim, P->obj->elements, P->obj->rowDim*P->obj->colDim*sizeof (double));
		
	return Q;
}
