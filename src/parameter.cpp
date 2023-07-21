#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameter.h"
#include "mystring.h"

static Parameter_t *parameter;

static void Parameter_load () {
	char 	fn[128] = "./Configuration";
	FILE*	fp = NULL;
	char 	buf[1024];
	char*	p = NULL;
	
	if ((fp = fopen (fn, "r")) == NULL) {
		fprintf (stderr, "to open the configuration file %s failed\n", fn);
		exit (0);
	}

	while (fgets (buf, sizeof (buf), fp)) {
		if (strlen(buf) < 2)
			continue;

		p = strtok (buf, " \t\n");
		if (p && strcmp (p, "ALGORITHM") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				strcpy (parameter->algorithm, p);
		} else  if (p && strcmp (p, "PROBLEM") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				strcpy (parameter->problem, p);
		} else 	if (p && strcmp (p, "NUMOBJ") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%d", &parameter->numObj);
		} else 	if (p && strcmp (p, "NUMVAR") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%d", &parameter->numVar);
		} else 	if (p && strcmp (p, "POPSIZE") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%d", &parameter->popSize);
		} else 	if (p && strcmp (p, "LIFETIME") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%lld", &parameter->lifetime);
		} else 	if (p && strcmp (p, "ID_CX") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%d", &parameter->id_cx);
		} else 	if (p && strcmp (p, "ID_MU") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%d", &parameter->id_mu);
		} else	if (p && strcmp (p, "RUN") == 0) {
			p = strtok (NULL, " \t\n");
			if (p)
				sscanf (p, "%d", &parameter->run);
		}
	}

	// close file 
	fclose (fp);
}



void Parameter_load (int argc, char **argv) {
	// 1. allocate memory for parameter 
	if ((parameter = (Parameter_t *)malloc (sizeof (Parameter_t))) == NULL) {
		fprintf (stderr, "%s:%d:allocate memory for parameter faild\n", __FILE__, __LINE__);
		exit (0);
	}

	// 2. get parameter from configuration file  
	Parameter_load ();

	// 3. get parameter from command line.
	if (argc <= 1) {

	} else if (argc < 8) {
		fprintf (stderr, "ERROR[%s:%d]: please use ./bin/main ALG PRO OBJ VAR Popsize FES RUN\n", __FILE__, __LINE__);
		exit (0);
	} else {
		// Algorithm
           	strcpy(parameter->algorithm, toUpper (argv[1]));

		// problem
                strcpy (parameter->problem, toUpper (argv[2]));

		// obj
               	parameter->numObj = atoi (argv[3]);

		// var
                parameter->numVar = atoi (argv[4]);

		// popSize 
               	parameter->popSize = atoi (argv[5]);

		// lifetime		
                parameter->lifetime = atoll (argv[6]);

		// run
                parameter->run = atoi (argv[7]);
	}

}

void Parameter_print (){
	if ( NULL == parameter) {
		fprintf (stderr, "ERROE:%s:%d: The parameter objective is NULL\n", __FILE__, __LINE__);
		exit (0);
	}
	printf ("Algorithm:\t%s\n", parameter->algorithm);
	printf ("Problem:\t%s\n", parameter->problem);
	printf ("numObj:\t\t%d\n", parameter->numObj);
	printf ("numVAR:\t\t%d\n", parameter->numVar);
	printf ("popSize:\t%d\n", parameter->popSize);
	printf ("lifetime:\t%lld\n", parameter->lifetime);
	printf ("id_cx:\t\t%d\n", parameter->id_cx);
	printf ("id_mx:\t\t%d\n", parameter->id_mu);
	printf ("run:\t\t%d\n", parameter->run);
}

Parameter_t *Parameter_get () {
	if ( NULL == parameter) {
		fprintf (stderr, "ERROE:%s:%d: The parameter objective is NULL\n", __FILE__, __LINE__);
		exit (0);
	}
	return parameter;
}
