#ifndef _HV_H
#define _HV_H

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "mystring.h"
#include "matrix.h"

/* 
   Top level of HV, search directory ./ouput/ ,
   calculate the HV and write it into a file.
*/
void HV ();

/*
   Get bound of matrix and write it into responding file.
*/
void hv_setBound ();

/*
Include and exclude algorithm
*/
double hv_iea (Matrix_t *M);

/*
The Lebesgue Algorithm calculates the hypervolume, Based on 
the paper: The measure of Pareto Optima: Application to Multiobjective
Metaheuristics
*/
double hv_leb (Matrix_t *M);


/*
 The algorithm HSO, based on the paper:
'A Faster Algorithm for Calculating HyperVolume'.
*/
double hv_hso (Matrix_t *M);

/*
 Based on 'On the Complexcity of Computing the Hypervolume Indicator'
*/
double hv_3d (Matrix_t *M);

/*
 According to 'A fast way of calculating exact hypervolumes' 
*/
double hv_wfg (Matrix_t *M);

/*
*/
double hv_2d (Matrix_t *M);

/*
 It has to be satisfied that M is a normilized front (non-dominated front).
*/
double hv (Matrix_t *M);

#endif
