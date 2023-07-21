#include "igd.h"
#include "problem.h"
#include "mystring.h"
#include "matrix.h"
#include "algebra.h"

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


void IGD () {
	struct 	dirent*	de = NULL;
	DIR*		dr = NULL;
	char		curfile[512];
	FILE*		fp = NULL;
	char**		items = NULL;
	char*		ptr = NULL;
	Matrix_t*	M = NULL ;
	double 		value;
	int 		i;

	// search the directory "./output"
	if ((dr = opendir ("./output")) == NULL) {
		fprintf (stderr, "Can't open the directory ./output \n");	
		exit (EXIT_FAILURE);
	}
	while ((de = readdir (dr)) != NULL) {
		if (strstr (de->d_name, "_front_") == NULL)
			continue;
		sprintf (curfile,"./output/%s",de->d_name);
		M = Matrix_read (curfile);

		sprintf (curfile, "%s", de->d_name);
		for (i=0; curfile[i] != '\0'; i++) {
			if (curfile[i] > '0' && curfile[i] <= '9') {
				curfile[i+1] = '\0';
				break;
			}
		}
		value = igd (M, curfile);

		// write the value of igd to a file
		sprintf (curfile,"./output/%s",de->d_name);
		items = str2item (curfile);
		ptr = strdup ("igd");	
		free (items[2]);
		items[2] = ptr;
		ptr = item2str (items);
		freeItem (items);

		if ((fp = fopen (ptr, "w")) != NULL) {
			fprintf (fp, "%lf\n", value);
			fclose (fp);
		} else {
			fprintf (stderr, "Can't open file %s\n", ptr);
			exit (EXIT_FAILURE);
		}
	}
}

static Matrix_t *igd_S = NULL;
double igd (Matrix_t *front, char *title) {
	double 		t;
	char		fn[1024];

	if (NULL == igd_S ) { 	
		// 1. generate sample points 
		igd_S = Problem_sample (title, front->colDim);

		// 2.
		if (NULL == igd_S)
			return 0;

		// 3. create the name of the PF file */
		sprintf (fn, "/tmp/PF_%s_%02d_%p", title, front->colDim, &t);

		/* 4. print the samples to file */
		Matrix_print (igd_S, fn);	
	} 

	t = distance_m2m (igd_S, front);
	return t;
}
