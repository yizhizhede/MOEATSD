#ifndef _TRANSITION_H
#define _TRANSITION_H

void b_poly (double *y, double alpha);
void b_flat (double *y, double A, double B, double C);
void b_param (double *y, double yy, double A, double B, double C);
void s_linear (double *y, double A);
void s_decept (double *y, double A, double B, double C);
void s_multi (double *y, double A, double B, double C);
double r_sum (double *y, double *w, int n);
double r_nonsep (double *y, int A, int n);

#endif
