#include "hv.h"
#include "dominate.h"

/*
  The algorithm is derived from paper "The Measure of Pareto Optima: Application to Multiobjective Metaheuristics"
  Suppose that the problem resolved is minization one, the desired result have been
  normalized. the reference point is set to (1.1, 1.1, ...).
*/

static double getBoundValue (Matrix_t *list, int i);
static void   spawnVector (Matrix_t *list, int i, double bi, Matrix_t **SpawnData);
static void   deleteP1 (Matrix_t *list);
static int    ndFilter (Matrix_t *list, Matrix_t *SpawnData);

double hv_leb (Matrix_t *M) 
{
	Matrix_t *list = Matrix_dup (M);
	Matrix_t *SpawnData = NULL;
	double 	LebMeasure = 0.0, lopOffVol = 0.0, bi = 0.0, lastVol = 0.0;
	int	newSize = list->rowDim, i;

	while (newSize > 1) {
		lopOffVol = 1.0;
		for (i=0; i < list->colDim; i++) {
			bi = getBoundValue (list, i);	
			spawnVector (list, i, bi, &SpawnData);
			lopOffVol *= (bi - list->elements[i]);
		}
		LebMeasure += lopOffVol;
		deleteP1 (list);
		newSize = ndFilter (list, SpawnData);
	}
	lastVol = 1.0;
	for (i=0; i < list->colDim; i++)
		lastVol *= (1.1 - list->elements[i]);
	LebMeasure += lastVol;
	return LebMeasure;
}

static double getBoundValue (Matrix_t *list, int i) {
	int row;
	double bi = 1.1, p1 = list->elements[i], t;

	for (row = list->rowDim -1; row > 0; row--) {
		t = list->elements[row * list->colDim + i];
		if ( t > p1 && t < bi) {
			bi = t;
		}
	}
	return bi;
}

static void spawnVector (Matrix_t *list, int i, double bi, Matrix_t **SpawnData) {
	int subset[] = {0};
	Matrix_t *subM = NULL;

	if ((*SpawnData) == NULL || (*SpawnData)->rowDim == 0) {
		(*SpawnData) = Matrix_sub (list, subset, 1);
		(*SpawnData)->elements[i] = bi;
	} else {
		subM = Matrix_sub (list, subset, 1);	
		subM->elements[i] = bi;
		Matrix_cat (SpawnData, subM);
		Matrix_free (&subM);
	}
}

static void deleteP1 (Matrix_t *list) {
	Matrix_t *dup = Matrix_dup (list);
	size_t size = (dup->rowDim -1) * dup->colDim * sizeof (double);

	memcpy (list->elements, dup->elements + dup->colDim, size);
	Matrix_free (&dup);
	list->rowDim--;
}

static int ndFilter (Matrix_t *list, Matrix_t *SpawnData) {
	int *flag = isDominated (SpawnData, list);	
	int i, j, rowDim = SpawnData->rowDim, colDim = SpawnData->colDim;
	int subSet[100], n=0;
	Matrix_t *subM = NULL;

	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			if (SpawnData->elements[i * SpawnData->colDim + j] - 1.1 > -1.0e-20)	{
				flag[i] = 1;
				break;
			}
		}
	}

	for (i=0, n=0; i<rowDim; i++) if (flag[i] == 0) {
		subSet[n++] = i;
	}
	free (flag);

	subM = Matrix_sub (SpawnData, subSet, n);
	Matrix_free (&SpawnData);

	// append Spawns to list
//	Matrix_cat (&list, subM);
//	Matrix_free (&subM);

	// Or prepend Spawns to list
	Matrix_cat (&subM, list);
	list->rowDim = subM->rowDim;
	free (list->elements);
	list->elements = subM->elements;
	free (subM);

	return list->rowDim;
}
