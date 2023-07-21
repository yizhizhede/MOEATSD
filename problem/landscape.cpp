#include "landscape.h"
#include "math.h"

#define PI 3.1415926535897932384626433832795028841971

double Sphere (double *x, int n) {
	double 	sum = 0;
	int 	i;

	for (i=0; i<n; i++) {
		sum += x[i] * x[i];
	}
	return sum;
}

double Schwefel (double *x, int n) {
	double 	max = 0;
	double 	t;
	int 	i;

	for (i=0; i<n; i++) {
		t = x[i];
		if (t < 0) 
			t = 0 - t;
		if (t  > max)
			max = t;
	}
	return max;
}

double Rosenbrock (double *x, int n) {
	double 	sum = 0;
	int 	i;

	for (i=0; i<n-1; i++) {
		sum += (100*(x[i]*x[i] - x[i+1])*(x[i]*x[i] - x[i+1]) + (x[i]-1)*(x[i]-1));
	}
	if (1 == n) {
		sum = (x[0]-1)*(x[0]-1);
	}
	return sum;
}

double Rastrigin (double *x, int n) {
	double 	sum = 0;
	int 	i;

	for (i=0; i<n; i++) {
		sum += (x[i]*x[i] - 10.0*cos(2*PI*x[i]) + 10.0);
	}
	return sum;
}

double Griewank (double *x, int n) {
	double 	sum = 0;
	double 	mul = 1.0;
	int 	i;

	for (i=0; i<n; i++) {
		sum += (x[i]*x[i]);
		mul *= cos (x[i]/sqrt (i+1));
	}
	sum = sum / 4000 - mul + 1.0;
	return sum;
}

double Ackley (double *x, int n) {
	double 	result;
	double 	sum1=0, sum2=0;
	int 	i;

	for (i=0; i<n; i++) {
		sum1 += (x[i] * x[i]);
		sum2 += (cos(2*PI*x[i]));
	}
	
	result = -20*exp(-0.2*sqrt(sum1/n)) - exp(sum2/n) + 20 + exp(1);
		
	return result;
}

// the landscape function of DTLZ
double dtlz_g1 (double *x, int n) {
	double t=0;
	int i;

	for (i=0, t=0.0; i<n; i++) {
		t += ((x[i]-0.5)*(x[i]-0.5) - cos(20*PI*(x[i]-0.5)));
	}
	return 100*(n+t);
}

double dtlz_g2 (double *x, int n) {
	double t=0;
	int i;	
	
	for (i=0; i<n; i++) {
		t += (x[i]-0.5)*(x[i]-0.5);
	}
	return t;
}

double dtlz_g6 (double *x, int n) {
	double t=0;
	int i;
	
	for (i=0; i<n; i++) {
		t += pow (x[i], 0.1);
	}
	return t;
}
 double dtlz_g7 (double *x, int n) {
	double t=0.0;
	int i;

	for (i=0; i<n; i++) {
		t += x[i];
	}
	return 1.0+9*t/n;
}
