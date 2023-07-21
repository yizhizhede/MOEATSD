#include "myrandom.h" 

int tournament(int n) {
	double c1, c2;

	c1 = n * randu ();
	c2 = n * randu ();

	if (c1 < c2)
		return (int)c1;
	else
		return (int)c2;
}
