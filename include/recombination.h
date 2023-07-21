#ifndef _RECOMBINATION_H
#define _RECOMBINATION_H

// 
void SBX_setPara (double index_mu, double index_cx);

// 
void   sbx_set_mutation (double index_mu);
double sbx_get_mutation ();

void   sbx_set_crossover (double index_cx);
double sbx_get_crossover ();

// 
double SBX_getNfitness ();

// reproduce : crossover + mutation
void SBX_reproduce (double* parent1, double* parent2, double* child1, double *f);

// mutation
void SBX_mutation (double* child1, double *f);


// crossover
void realbinarycrossover(double *parent1, double *parent2, double *child1, double *child2,
		                           double rate, int len, double *low, double *upp);
// mutation
void realmutation(double* ind, double rate, int len, double *low, double *upp);

#endif
