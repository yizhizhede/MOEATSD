#ifndef _MYRANDOM_H 
#define _MYRANDOM_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


/* generate distribute function sample  */
double randnormal(double mu,double sigma);
double randn();

double randlevy(double mu,double scale);
double randl(double mu);

double randcauchy(double mu,double c);
double randc();

//
double randu();

//
int randmatrix(double *matrix,int nrows,int ncols);

// exponetial distribution
double randexp (double lambda);

#endif 
