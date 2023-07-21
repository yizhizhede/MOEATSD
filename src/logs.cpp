#include "logs.h"
#include <string.h>

// #define LOG_START

void log_string (char *fn, char *str) {
#ifdef LOG_START
	FILE*	fp = NULL;
	char 	buff[1024];

	if (fn == NULL || strlen (fn) == 0)	
		return;

	sprintf (buff, "./log/%s", fn);
	fp = fopen (buff, "a");
	fprintf (fp, "%s\n", str);
	fclose (fp);
#endif
}
void log_matrix (char *fn, Matrix_t *M) {
#ifdef LOG_START
	char 	buff[1024];
	FILE*	fp = NULL;
	int 	i, j;

	if (fn == NULL || strlen (fn) == 0)	
		return;
	sprintf (buff, "./log/%s", fn);
	fp = fopen (buff, "a");
	for (i=0; i<M->rowDim; i++) {
		for (j=0; j<M->colDim; j++) {
			fprintf (fp, "%f ", M->elements[i*M->colDim+j]);
		}
		fprintf (fp, "\n");
	}
	fprintf (fp, "\n");
	fclose (fp);
#endif
}

void log_vector (char *fn, int *buf, int len) {
#ifdef LOG_START
	FILE*	fp = NULL;
	char 	buff[1024];
	int	i;

	if (fn == NULL || strlen (fn) == 0)	
		return;

	sprintf (buff, "./log/%s", fn);
	fp = fopen (buff, "a");
	for (i=0; i<len; i++) {
		fprintf (fp, "%d ", buf[i]);
	}	
	fprintf (fp, "\n");
	fclose (fp);
#endif
}

void log_vector (char *fn, double* buf, int len) {
#ifdef LOG_START
	FILE*	fp = NULL;
	char 	buff[1024];
	int	i;

	if (fn == NULL || strlen (fn) == 0)	
		return;

	sprintf (buff, "./log/%s", fn);
	fp = fopen (buff, "a");
	for (i=0; i<len; i++) {
		fprintf (fp, "%f ", buf[i]);
	}	
	fprintf (fp, "\n");
	fclose (fp);
#endif
}

