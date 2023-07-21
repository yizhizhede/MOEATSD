#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "matrix.h"

typedef struct Problem_tag {
	char 	title[32];	/* 1. The title of the problem */

	int	numObj;		/* 2. The number of decision variables */
	int 	numVar;		/* 3. The number of objective functions */
	
	double 	*lowBound;	/* 4. The lower boundary of decision variables */
	double 	*uppBound;	/* 5, The upper boundary of decision variables */

				/* 6. pointer of fitnenss function */ 
	void 	(*evaluate) (double *var, int numVar, double *obj, int numObj); 	

	long long fitness;	/* 7. The counter of fitness function */ 
	long long lifetime;	/* 8. The terminal condition */

	double	*ideal_point;	/* 9. The ideal point in objective space */
} Problem_t;

Problem_t* 	Problem_new ();
Problem_t* 	Problem_get ();
Matrix_t*  	Problem_sample (char *title, int numObj);
void       	Problem_evaluate (double *var, int numVar, double *obj, int numObj);
void       	Problem_progressbar ();
int        	Problem_isTick ();
int 	   	Problem_isEnd ();
long long  	Problem_getFitness ();
long long  	Problem_getLifetime ();
double*   	Problem_getLowerBound ();
double*   	Problem_getUpperBound ();
double*   	Problem_getIdealPoint ();
void	   	Problem_xToOne (Matrix_t* X);
void    	Problem_xFromOne (Matrix_t* X);
void	   	Problem_xToOne (double* X, int numVar);
void    	Problem_xFromOne (double*  X, int numVar);

#endif
