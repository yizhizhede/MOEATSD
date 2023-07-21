#ifndef _MYSTRING_H
#define _MYSTRING_H

char**	str2item (char*  str);
char*	item2str (char** list);
void 	showItem (char** list);
void 	freeItem (char** list);
char*	toUpper  (char* str);
char*	strrep   (char* src, char* oldStr, char *newStr);

#endif
