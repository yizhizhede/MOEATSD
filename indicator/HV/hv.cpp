#include "hv.h"
#include "link.h"
#include "dominate.h"

//#define HV_LEB 	0
//#define HV_IEA 	1
//#define HV_HSO 	2
//#define HV_3D 	3
//#define HV_WFG  	4

/*
 It has to be satisfied that M is a normilized front (non-dominated front). 
 otherwise it is not promised to get a correct answer.
 And the Matrix M is normalized by the maxmum and minmum, the reference point 
 is set to (1.1, 1.1, ...) in default.
*/
double hv (Matrix_t *M) {
	// 1. 
	if (M == NULL || M->elements == NULL || M->rowDim < 1 || M->colDim < 1) {
		return 0;
	}

	// 2. 
	switch (M->colDim) {
		case 2:
			return hv_2d (M);	
	//	case 3:
	//		return hv_3d (M);	
		default:
			return hv_wfg (M);
	}
	return 0;
}
 


void HV ()
{
	struct 	dirent *de = NULL;
	DIR 	*dr = NULL;
	char	curfile[512];
	FILE	*fp = NULL;
	char 	**items = NULL;
	char 	*ptr = NULL;
	Matrix_t *M = NULL, *bound = NULL, *normM = NULL, *F1 = NULL;
	List_t *list=NULL;
	int 	*arr=NULL;
	double 	hv = -1.0;

	// setBound
	hv_setBound ();

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

		items = str2item (curfile);
		sprintf (curfile, "%s_bound", items[0]);
		bound = Matrix_read (curfile);

		normM = Matrix_norm (M, bound);	
		Matrix_free (&M);
		Matrix_free (&bound);

		list = ndSort (normM);	
		arr = Link2Array (list->list_head);
		F1 = Matrix_sub (normM, arr+1, arr[0]);
		Matrix_free (&normM);
		free (arr);
		List_free (&list);
		

		// specific hv function
#ifdef HV_LEB
		printf ("hv_leb: %lf\n", hv_leb (F1));
#endif
#ifdef HV_IEA
		printf ("hv_iea: %lf\n", hv_iea (F1));
#endif
#ifdef HV_HSO
		printf ("hv_hso: %lf\n", hv_hso (F1));
#endif
#ifdef HV_3D
		if (F1->colDim ==3)
			printf ("hv_3d: %lf\n", hv_3d (F1));
#endif
#ifdef HV_WFG
		printf ("hv_wfg: %lf\n", hv = hv_wfg (F1));
#endif

		// printf ("%s\n", de->d_name);
		hv = hv_wfg (F1);
		Matrix_free (&F1);

		// write hv to a file
		ptr = strdup ("hv");	
		free (items[2]);
		items[2] = ptr;
		ptr = item2str (items);
		freeItem (items);

		if ((fp = fopen (ptr, "w")) != NULL) {
			fprintf (fp, "%lf\n", hv);
			fclose (fp);
		} else {
			fprintf (stderr, "Can't open file %s\n", ptr);
			exit (EXIT_FAILURE);
		}
	}
}
