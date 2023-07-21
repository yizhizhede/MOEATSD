#include "gd.h"
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


void GD () {
	struct 	dirent *de = NULL;
	DIR 	*dr = NULL;
	char	curfile[512];
	FILE	*fp = NULL;
	char 	**items = NULL;
	char 	*ptr = NULL;
	Matrix_t *M = NULL ;
	double 	value;
	int i;

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
		value = gd (M, curfile);

		// write the value of gd to a file
		sprintf (curfile,"./output/%s",de->d_name);
		items = str2item (curfile);
		ptr = strdup ("gd");	
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

double gd (Matrix_t *front, char *title) {
	double t;
	Matrix_t *S = Problem_sample (title, front->colDim);
	if (S==NULL)
		return 0;
	t = distance_m2m (front, S);
	Matrix_free (&S);
	return t;
}
