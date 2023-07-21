
/* Box_mullter algorithm generate Gaussian distribute sample
 * https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform */
#include "myrandom.h"

#define PI 3.14159265358979323846

double 
randn()
{
	double u1,u2;

	const double epsilon=1.0/pow(10.0,15);
	const double two_pi=2.0*3.14159265358979323846;
	
	static double z0,z1;
	static int generate=1;
	generate = !generate;

	if (!generate) 
		return z1;

	do {
		u1 = rand()*(1.0/ RAND_MAX);
		u2 = rand()*(1.0/ RAND_MAX);
	}
	while (u1 <= epsilon);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

	return z0;
}

double 
randnormal(double mu,double sigma)
{
	return mu+sigma*randn();
}

double 
randc()
{
	double u;
	do {
		u=(double)rand()/RAND_MAX;
	} while(u==0 || u==1);
	u=PI*(u-0.5);
	return tan(u);
}
double
randcauchy(double mu,double c)
{
	return c*randc()+mu;
}

/* Levy distribution 
 * http://markread.info/2016/08/code-to-generate-a-levy-distribution/  
 */
/* Samples a levy distribution wherein the power law decay can be 
 * adjusted between 1/x and 1/x^3.
 * Note that this sampling method can return negative value. Values 
 * are symmetrical around zero.
 * param mu must lie between 1 and 3. Correspods to 1/x and 1/x^3*/
double
sample(double mu)
{
	double u,v,t,s;
	double alpha=mu-1;
		
	if ( (mu<=1)||(mu>3))
	return 0;

	do {
		u=(double)rand() / RAND_MAX;	
		v=(double)rand() / RAND_MAX;
	} while(!v);				// discard 0
	u=PI*(u-0.5);				// -PI/2 < u < PI/2	
	v=-log(v);				

	//general case
	t =  sin(alpha*u)/pow(cos(u),1/alpha);
	s =  pow(cos((1-alpha)*u)/v,(1-alpha)/alpha);
	return t*s;
}

/* Same as,but ensure all values are positive. Negative value are simple
 * negated, as the levy distribution represented is symmetrical 
 * around zero .*/
double 
randlevy(double mu,double scale)
{
	double l=scale*sample(mu);
	if (l<0.0) 
		return -1.0 * l;
	return l;
}

/* Default value case,scale=1*/
double 
randl(double mu)
{
	return randlevy(mu,1.0);
}

double 
randu()
{
	return 1.0 * rand() / (RAND_MAX + 1.0);
}

int 
randmatrix(double *matrix,int nrows,int ncols)
{
	int i;
	
	for (i=nrows*ncols-1; i>=0; i--)
		matrix[i] = randu();

	return 0;
}

double randexp (double lambda) {
	double U = randu ();

	if (U == 0)
		return 1.0 / lambda;
	else
		return (-1.0 / lambda) * log(U);
}
