#include "interactive.h"
#include "problem.h"
#include "myrandom.h"
#include <string.h>
#include <float.h>

int isInteractive (int i, int j) {
	return isInteractive (i, j, 6);
}

int isInteractive (int i, int j, int NIA) {
	Problem_t*	problem = Problem_get ();
	int		numVar  = problem->numVar;
	int		numObj  = problem->numObj;
	double*		lowBound = problem->lowBound;
	double*		uppBound = problem->uppBound;
	double 		x[4*numVar];
	double		y[4*numObj];
	int		k, m;
	double		a2, b2, delta1, delta2;

	for (k=0; k<NIA; k++) {
		for (m=0; m<numVar; m++) {
			x[m] = lowBound[m] + randu() * (uppBound[m] - lowBound[m]);
		}
		memcpy (x+1*numVar, x, numVar*sizeof (double));
		memcpy (x+2*numVar, x, numVar*sizeof (double));
		memcpy (x+3*numVar, x, numVar*sizeof (double));
		a2 = lowBound[i] + randu() * (uppBound[i] - lowBound[i]);
		b2 = lowBound[j] + randu() * (uppBound[j] - lowBound[j]);
		x[1*numVar+i] = a2;
		x[2*numVar+j] = b2;
		x[3*numVar+i] = a2;
		x[3*numVar+j] = b2;
		Problem_evaluate (x+0*numVar, numVar, y+0*numObj, numObj);
		Problem_evaluate (x+1*numVar, numVar, y+1*numObj, numObj);
		Problem_evaluate (x+2*numVar, numVar, y+2*numObj, numObj);
		Problem_evaluate (x+3*numVar, numVar, y+3*numObj, numObj);
		for (m=0; m<numObj; m++) {
			delta1 = y[1*numObj+m] - y[0*numObj+m];
			delta2 = y[3*numObj+m] - y[2*numObj+m];
			
			printf ("i, j = %d, %d\n", i, j);
			printf ("x0 = %f %f, %f %f\n", x[0*numVar+i], x[0*numVar+j], y[0*numObj+0], y[0*numObj+1]);
			printf ("x1 = %f %f, %f %f\n", x[1*numVar+i], x[1*numVar+j], y[1*numObj+0], y[1*numObj+1]);
			printf ("x2 = %f %f, %f %f\n", x[2*numVar+i], x[2*numVar+j], y[2*numObj+0], y[2*numObj+1]);
			printf ("x3 = %f %f, %f %f\n", x[3*numVar+i], x[3*numVar+j], y[3*numObj+0], y[3*numObj+1]);

			printf ("delta=%.16f, k=%d, m=%d\n", fabs(delta1 - delta2), k, m);
			// 1. base on DG and DG2
			if (fabs(delta1 - delta2) > 2*DBL_EPSILON) {
				return 1;
			}

			// 2. base on MOEA/DVA
			// if (delta1*delta2 < 0) {
			//	return 1;
			// }
		}
	}
	return 0;
}
