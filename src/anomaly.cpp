#include "anomaly.h"
#include <math.h>
#include <float.h>


int anomaly (double* A, int n) {
	int 	i, k; 
	double	miu = 0, delta = 0;
	double  Z[n+10];
	double 	z;

	for (i=0; i<n; i++) {
		miu += A[i];
	}
	miu /= n;

	for (i=0; i<n; i++) {
		delta += (A[i] - miu)*(A[i] - miu);
	}
	delta = sqrt (delta / n);
		
	if (delta < DBL_EPSILON)
		return -1;
	for (i=0; i<n; i++) {
		Z[i] = (A[i] - miu) / delta;
	}
	for (i=1, k=0, z=Z[0]; i<n; i++) {
		if (fabs(z) < fabs(Z[i])) {
			k = i;
			z = Z[i];
		}
	}
	if (z > 2.5 || z < -2.5) {
		return k;
	}
	return -1;
}
