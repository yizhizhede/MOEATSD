#include "kmeans.h"
#include "algebra.h"

#include <string.h>
#include <stdlib.h>

int *kmeans(double *Pop, int nrows, int ncols)
{
	int 	i,j;
	double 	c[4 * ncols];
	int 	c1,c2;
	int 	t;
	int 	cluster[nrows];
	double 	temp;
	double 	tt;
	double 	d1,d2;
	int 	n=2;
	int *	flag = (int *)malloc (nrows * sizeof (int));

	// for terminal
	int terminal = 0;

	// initial center pointer
	t = 1.0 * nrows * rand() / RAND_MAX;
	memcpy(c, Pop + t*ncols, ncols*sizeof(double));

	for (i=0,temp=0.0; i<nrows; i++) {
		tt = distance_p2p(Pop+i*ncols, c, ncols); 
		if (tt > temp) {
			t = i;
			temp = tt;
		}
	}
	memcpy(c+ncols, Pop+t*ncols, ncols*sizeof(double));

	//
	do {
		n = (n+2) % 4;
		c1=c2=0;
		memset(c+((n+2)%4)*ncols, 0, 2*ncols*sizeof(double));

		for (i=0; i<nrows; i++) {
			d1 = distance_p2p(Pop+i*ncols,c+n*ncols,ncols);
			d2 = distance_p2p(Pop+i*ncols,c+(n+1)*ncols,ncols);
			if (d1 <= d2) {
				cluster[i] = 1;
				for (j=0; j<ncols; j++) {
					c[((n+2)%4)*ncols+j] += Pop[i*ncols+j];
				}
				c1++;
			}
			else {
				cluster[i] = 2;
				for (j=0; j<ncols; j++) {
					c[((n+3)%4)*ncols+j] += Pop[i*ncols +j];
				}
				c2++;
			}
		}
		
		for (j=0; j<ncols; j++) {
			c[((n+2)%4)*ncols + j ] /= c1;
			c[((n+3))%4*ncols + j ] /= c2;
		}

		// 
		for (i=0; i<2*ncols; i++) if (c[i]!=c[i+2*ncols]) break;
		if (i==2*ncols) 
			t = 0;
		else 
			t =1;
		
		// for terminal
		terminal++;
		if (terminal>20)
			break;
	} while(t);
	
	//
	d1 = sum(c,ncols);
	d2 = sum(c+ncols,ncols);

	if (d1 < d2) {
		for (i=0; i<nrows; i++) {
			if (cluster[i]==1) {
				flag[i] = 0;
			} else {
				flag[i] = 1;
			}
		}
	} else {
		for (i=0; i<nrows; i++) {
			if (cluster[i] ==2) {
				flag[i] = 0;
			} else {
				flag[i] =1;
			}
		}
	}

	return flag;
}
