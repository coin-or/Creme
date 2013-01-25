/*
 * Name:    readbz.c
 * Author:  Pietro Belotti
 * Purpose: check that bzgetXXX functions work
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 */

#include <bzlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lpio.h"
#include "misc.h"

/* #define MAX_STR 50000 */

/**********************************/

int main (int argc, char **argv) {

  int i;

  for (i=1; i<argc; i++) {

    int status;
    double x;
    /*char c;*/

    BZFILE *bzf;
    FILE *f;

    fprintf (stderr, "reading %s\n", argv [i]);

    f = fopen (argv [i], "r");

    bzf = BZ2_bzReadOpen( &status, f, 1, 0, NULL, 0);

    do {
      status = bzgetdbl (bzf, &x);
      printf ("[%.10f]\n", x);
    } while  (!status);

    /*
    do {
      status = bzgetchar (bzf, &c);
      if (!status) printf ("%c", c);
    } while  (!status);
    */

    BZ2_bzReadClose (&status, bzf);
    fclose (f);
  }

  return 0;
}
