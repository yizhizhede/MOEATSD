#include "hv.h"

#include <stdlib.h>
#include <string.h>

double hv_2d (Matrix_t *M) {
	int 	*index = sort (M, 0);
	double	volume = 0; 
	int i;
	
	volume = (1.1 - M->elements[2*index[0]])*(1.1-M->elements[2*index[0]+1]);
	for (i=1; i<M->rowDim; i++) {
		volume += (1.1 - M->elements[2*index[i]])*(M->elements[2*index[i-1]+1]-M->elements[2*index[i]+1]);
	}

	free (index);
	return volume;
}
