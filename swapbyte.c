#include "stdio.h"

void swapbyte4_(field, nx)
int *field;
int *nx;
{
    /* Local variables */
    static int j;
    int t;
    char *p,*q;

    /* Function Body */
    for (j = 0; j < *nx; ++j) {
	p = (char *) (field + j);
	p += 3;
	q = (char *) &t;
	*q++ = *p--;
	*q++ = *p--;
	*q++ = *p--;
	*q++ = *p--;
	field[j] = t;
    }
} /* swapbyte_ */

