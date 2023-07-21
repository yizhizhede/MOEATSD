#include "transition.h"
#include <math.h>
#include <string.h>
#include <float.h>

#define PI 3.14159265358979323846264338327950288419716939937510

#define MAX(a,b) ((a>b)?(a):(b))
#define MIN(a,b) ((a<b)?(a):(b))


void b_poly (double *y, double alpha) {
	if (*y < DBL_EPSILON) {
		*y = 0;
	} else {
		*y = pow (*y, alpha);
	}
}

void b_flat (double *y, double A, double B, double C) {
	*y = A + MIN(0, floor (*y - B))*A*(B-*y)/B - MIN(0,floor (C - *y))*(1-A)*(*y-C)/(1-C);
}

void b_param (double *y, double u, double A, double B, double C) {
	double v = A-(1.0-2.0*u)*fabs(floor (0.5-u)+A);

	if (*y < DBL_EPSILON) {
		*y = 0;
	} else {
		*y = pow(*y, B+(C-B)*v);
	}
}

void s_linear (double *y, double A) {
	*y = fabs(*y-A)/fabs(floor(A-*y)+A);
}

void s_decept (double *y, double A, double B, double C) {
	*y = 1.0+(fabs(*y-A)-B)*
		(floor(*y-A+B)*(1.0-C+(A-B)/B)/(A-B) + floor(A+B-*y)*(1.0-C+(1.0-A-B)/B)/(1.0-A-B) + 1.0/B);
} 

void s_multi (double *y, double A, double B, double C) {
	double t = fabs(*y-C)/(2.0*(floor(C-*y)+C));
	*y = (1.0+cos((4*A+2)*PI*(0.5-t))+4.0*B*t*t)/(B+2.0);
}

double r_sum (double *y, double *w, int n) {
	double 	t1=0, t2=0;
	int 	i;

	for (i=0; i<n; i++) {
		t1 += w[i]*y[i];
		t2 += w[i];
	}

	return t1/t2;
}

static int cmp (const void *a, const void *b) {	
	double 	t;

	t = *(double *)a - *(double *)b; 
	if (t > 0)
		return 1;
	else
		return 0;
}//升序

static double r_nonsep (double *y, int n) {
	double 	z[n+10];
	int 	j, A = n;
	double 	t1=0, t2=0;
	
	memcpy (z, y, n*sizeof (double));
	qsort (z, n, sizeof(double), cmp);

	for (j=1; j<n; j++) {
		t2 += 2*j*(n-j)*(z[j]-z[j-1]);
	}
	for (j=0; j<n; j++) {
		t1 += z[j];
	}
	t1 += t2;

	return t1/((1.0*n/A)*ceil(A/2.0)*(1.0+2*A-2.0*ceil(A/2.0)));
}

double r_nonsep (double *y, int A, int n) {
	int 	j, k;
	double 	t1=0, t2=0;

	if (A == n) {
		return	r_nonsep (y, n);
	} else {
		for (j=0; j<n; j++) {
			t1 += y[j];
			for (k=0, t2=0.0; k<=(A-2); k++) {
				t2 += fabs(y[j]-y[(j+1+k)%n]);
			}
			t1 += t2;
		}
		return t1/((1.0*n/A)*ceil(A/2.0)*(1.0+2*A-2.0*ceil(A/2.0)));
	}
}

