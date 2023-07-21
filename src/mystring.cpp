#include "mystring.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


char **str2item (char *str) {
	if (str == NULL || strlen (str) < 1) {
		fprintf (stderr, "string is NULL or empty\n");
		exit (EXIT_FAILURE);
	}

	char *tofree, *string, *tok;
	int i = 0;
	char **list = (char **)malloc (110 * sizeof (char *));

	tofree = string = strdup (str);
	
	i = 0;
	while ((tok = strsep (&string, "_")) != NULL && i < 100) {
		list[i++] = strdup (tok);
	}
	list[i] = NULL;
	free (tofree);

	return list;
}

char *item2str (char **list) {
	int 	len = 0, i;
	char*	str = NULL;
	char*	buf = NULL;
	
	for (i=0; i < 100 && list[i] != NULL; i++) {
		len += strlen (list[i]);	
	}

	len += 10;
	str = (char *)malloc (len * sizeof (char ));
	buf = (char *)malloc (len * sizeof (char ));
	strcpy (str, list[0]);

	for (i=1; i < 100 && list[i] != NULL; i++) {
		strcpy (buf, str);
		sprintf (str, "%s_%s", buf, list[i]);
	}
	free (buf);
	
	return str;
}

void showItem (char **list) {
	int 	i=0;
	
	for (i=0; i < 100 && list[i] != NULL; i++) {
		printf ("%s\n", list[i]);
	}
}

void freeItem (char **list) {
	int i=0;
	
	for (i=0; i < 100 && list[i] != NULL; i++) {
		free(list[i]);
	}
}

char *toUpper (char *str) {
	int 	i;
	char*	upp = (char *)malloc ((strlen (str)+1)*sizeof (char));

	if (str==NULL)
		return NULL;

	for (i=0; str[i] != '\0'; i++) {
		if (str[i] >= 'a' && str[i] <= 'z') {
			upp[i] =(char)(str[i] - 'a' + 'A');
		} else {
			upp[i] = str[i];
		}
	}

	upp[i] = '\0';
	return upp;
}


char* strrep (char *src, char* oldStr, char *newStr) {
	int 	len1 = strlen (src);
	int	len2 = strlen (oldStr);
	int	len3 = strlen (newStr);
	char*	buf  = NULL;
	int	i, j, k;

	if (src == NULL) 
		return NULL;
	if (oldStr == NULL || newStr == NULL) 
		return src;

	buf = (char *)malloc ((len1+len3+10)*sizeof (char));

	for (i=0; i<len1; i++) {
		for (j=0; j<len2; j++) {
			if (src[i+j] != oldStr[j]) {
				break;
			}
		}
		if (j >= len2) {
			break;
		}
	}
	if (i>=len1) 
		return src;

	strcpy (buf, src);
	for (j=0, k=i; j<len3; j++) {
		buf[k++] = newStr[j];
	}
	for (j=i+len2; j<len1; j++) {
		buf[k++] = src[j];
	}
	buf[k++] = '\0';

	return buf;	
}
