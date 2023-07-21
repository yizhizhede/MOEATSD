#include "recombination.h"
#include "problem.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

static double 	id_mu = 20;
static double 	id_cx = 20;

// create a random U \in [0, 1) 
static  double 	next_U ();


double SBX_getNfitness () {
	return Problem_getFitness ();
}

void SBX_setPara (double index_mu, double index_cx) {
	id_mu = index_mu;
	id_cx = index_cx;	
}

void sbx_set_mutation (double index_mu) {
	id_mu = index_mu;
}

double sbx_get_mutation () {
	return id_mu;
}

void sbx_set_crossover (double index_cx) {
	id_cx = index_cx;	
}

double sbx_get_crossover () {
	return id_cx;	
}

void SBX_reproduce (double* parent1, double* parent2, double* child1, double *f) {
	Problem_t *problem = Problem_get ();
	int len = problem->numVar;
	double child2[len];	

	realbinarycrossover(parent1, parent2, child1, child2, 1.0, len, problem->lowBound, problem->uppBound);
	realmutation(child1, 1.0/len, len, problem->lowBound, problem->uppBound);

	Problem_evaluate (child1, len, f, problem->numObj);
}

void SBX_mutation (double* child1, double *f) {
	Problem_t *problem = Problem_get ();
	int len = problem->numVar;

	realmutation(child1, 1.0/len, len, problem->lowBound, problem->uppBound);

	Problem_evaluate (child1, problem->numVar, f, problem->numObj);
}

void realmutation(double* ind, double rate, int len, double *low, double *upp)
{
	double 	rnd, delta1, delta2, mut_pow, deltaq;
	double 	y, yl, yu, alpha;
	double 	eta_m = id_mu;
	int 	i;

	for (i=0; i<len; i++) { 
		if((rnd = next_U()) <= rate) {
			y = ind[i];
			yl = low[i];
			yu = upp[i];
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);

			rnd = next_U(); 
			mut_pow = 1.0/(eta_m + 1.0);
			if (rnd <= 0.5) {
				alpha = 2.0*rnd + (1.0 - 2.0*rnd)*pow (1.0 - delta1, eta_m + 1.0);
				deltaq = pow (alpha, mut_pow) - 1.0;
			} else {
				alpha = 2.0*(1.0-rnd) + 2.0*(rnd - 0.5)*pow (1.0 - delta2, eta_m + 1.0);
				deltaq = 1.0 - pow (alpha, mut_pow);
			}
			y = y + deltaq*(yu - yl);

			if (y > yu)
				y = yu;
			if (y < yl)
				y = yl;
			ind[i] = y;
		}
	}
}


void realbinarycrossover(double *parent1, double *parent2,double *child1,double *child2, 
						double rate, int len, double *low,double *upp) 
{
	double 	rnd;
	double 	y1, y2, yl, yu;
	double 	c1, c2;
	double 	alpha, beta, betaq;
	double 	eta_c = id_cx;
	int 	i;

	if ((rnd = next_U ()) <= rate) {
		for (i=0; i<len; i++) {
			if ((rnd = next_U()) <= 0.5) {
				if (fabs(parent1[i] - parent2[i]) > DBL_EPSILON) {
					if (parent1[i] < parent2[i]) {
						y1 = parent1[i];
						y2 = parent2[i];
					} else {
						y1 = parent2[i];
						y2 = parent1[i];
					}
					yl = low[i];
					yu = upp[i];
					rnd = next_U ();

					// compute c1
					beta = 1.0 + 2.0*(y1-yl)/(y2 - y1);
					alpha = 2.0 - pow (beta, -(eta_c + 1.0));
					if (rnd <= (1.0/alpha)) {
						betaq = pow (rnd*alpha, 1.0/(eta_c + 1.0));
					} else {
						betaq = pow (1.0/(2.0 - rnd*alpha), 1.0/(eta_c + 1.0));
					}
					c1 = 0.5 * ((y1 + y2) - betaq*(y2-y1));

					// compute c2
					beta = 1.0 + 2.0*(yu - y2)/(y2 - y1);
					alpha = 2.0 - pow (beta, -(eta_c + 1.0));
					if (rnd <= (1.0/alpha)) {
						betaq = pow (rnd*alpha, 1.0/(eta_c + 1.0));
					} else {
						betaq = pow (1.0/(2.0 - rnd*alpha), 1.0/(eta_c + 1.0));
					}	
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));

					// make sure yl <= c1, c2 <= yu.
					if (c1 < yl)
						c1 = yl;
					if (c2 < yl)
						c2 = yl;
					if (c1 > yu)
						c1 = yu;
					if (c2 > yu)
						c2 = yu;

					// set the value of child into c1, c2.
					if ((rnd = next_U()) <= 0.5) {
						child1[i] = c2;
						child2[i] = c1;
					} else {
						child1[i] = c1;
						child2[i] = c2;
					}
				} else {
					child1[i]=parent1[i];
					child2[i]=parent2[i];
				}
			} else {
				child1[i]=parent1[i];
				child2[i]=parent2[i];
			}
		}
	} else {
		for (i=0; i<len; i++) {
			child1[i] = parent1[i];
			child2[i] = parent2[i];
		}
	}	
}

#define	Q_SBX_LEN 1000000
static 	int 	Q_sbx_index = Q_SBX_LEN;
static 	double 	Q_sbx[Q_SBX_LEN+10];
static  double 	next_U () {
	int 	i;
	double	U;

	if (Q_sbx_index >= Q_SBX_LEN) {
		Q_sbx_index = 0;
		srand (time (NULL));
		for (i=0; i<Q_SBX_LEN; i++) {
			Q_sbx[i] = 1.0 * rand () / (RAND_MAX + 1.0);
		}
	}

	U = Q_sbx[Q_sbx_index++];
	return U;
}
