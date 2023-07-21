#include "hv.h" 
#include "dominate.h"
#include "link.h"

static void hv_removeBound ()
{
	struct 	dirent *de = NULL;
	DIR	*dr = NULL;
	char 	path[512]="\0";	

	dr = opendir ("./output");
	if (dr == NULL) {
		fprintf (stderr, "Can't open ./output directory\n");
		exit (EXIT_FAILURE);
	}
	while ((de = readdir (dr)) != NULL) {
		if (strstr (de->d_name, "bound")) {
			sprintf (path, "./output/%s", de->d_name);
			if (remove (path) == -1) {
				fprintf (stderr, "Remove file %s failure\n", path);
			}
		}
	}
}

void hv_setBound ()
{
	struct 	dirent *de = NULL; 
	DIR 	*dr = NULL;
	struct	stat sb; 
	char 	curfile[512];
	char 	**items = NULL;
	Matrix_t *M = NULL, *bound = NULL, *T=NULL;
	List_t *list = NULL;
	int *arr = NULL;

	// remove the bound file created by previous procedure
	hv_removeBound ();

	// file = ./output/pro_alg_con_id
	dr = opendir ("./output");	
	if (dr == NULL) {
		fprintf (stderr, "Can't open current directory\n");
		exit (EXIT_FAILURE);
	}
		
	while ((de = readdir (dr)) != NULL) {
		if (strcmp (de->d_name,".") == 0 || strcmp (de->d_name, "..") == 0)
			continue;
		if ( strstr (de->d_name, "front") == NULL) 
			continue;
		strcpy (curfile, "./output/");
		strcat (curfile, de->d_name);
		if (stat (curfile, &sb) == -1) {
			perror ("stat");
			exit (EXIT_FAILURE);
		}
		if ((sb.st_mode & S_IFMT) == S_IFREG) {
			T = Matrix_read (curfile);
			list = ndSort (T);
			arr = Link2Array (list->list_head);
			M = Matrix_sub (T, arr+1, arr[0]);

			Matrix_free (&T);
			List_free (&list);
			free (arr);

			bound = Matrix_bound (M);	
		
			items = str2item (curfile);
			sprintf (curfile, "%s_bound", items[0]);
			freeItem (items);

			if (stat (curfile, &sb) == -1) { 	// file does not exist.
				Matrix_print (bound, curfile);
			} else { 				// file exists
				Matrix_free (&M);
				M = Matrix_read (curfile);
				Matrix_cat (&M, bound);
	
				Matrix_free (&bound);
				bound = Matrix_bound (M);
				Matrix_print (bound, curfile);
			}
		}
	}

	Matrix_free (&M);
	Matrix_free (&bound);

	closedir (dr);
}
