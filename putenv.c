#include <stdlib.h>
#include <string.h>
#include "f2c.h"

integer putenv_(string, string_len)
char *string;
ftnlen string_len;
{
    /* locals */
    integer i,ret_val;
    char *newstring;

    newstring = malloc( (string_len + 1) * sizeof(char) );
    newstring = strncpy(newstring,string,string_len);
    newstring[string_len] = '\0';
    i = string_len - 1;
    while ( newstring[i] == ' ' && i > 0 ) {
        newstring[i--] = '\0';
    }
    ret_val = putenv(newstring);
    free(newstring);
    return ret_val;
} /* putenv_ */

