/*
 * Name:    writelp.c
 * Author:  Pietro Belotti
 * Purpose: write feasible subsystem into a .lp file
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <bzlib.h>

#include "sparse.h"
#include "lpio.h"
#include "rtr.h"
#include "misc.h"

#include "config_creme.h"

/*
 * Write feasible subsystem into a .lp file
 */

int write_lp (sparseLP *lp,   /* LP data   */
	      char *mfs,      /* feasible subsystem */
	      char *filename  /* file name */
	      ) {

  FILE *f;

  register int i;

  if (strcmp (filename, "-")) f = fopen (filename, "w");
  else                        f = stdout;

  if ((f == NULL) || ferror (f)) {
    printf ("Unable to open file %s\n", filename);
    return -1;
  }

  fprintf (f, "\\ Feasible subsystem created by Crème %d.%d\n\n"
	   "minimize x0\n\n"
	   "subject to\n\n", 
	   CREME_VERSION_MAJOR, CREME_VERSION_MINOR);

  for (i=0; i<lp -> rk; ++i)
    if (mfs [i]) {
      int j;
      fprintf (f, "c%04d: ", i);
      for (j=0; j < lp -> il [i]; ++j)
	fprintf (f, "%+g x_%d ", lp -> ic [i][j], lp -> ip [i][j]);
      fprintf (f, ">= %e\n", lp -> rhs [i]);
    }

  fprintf (f, "\nbounds\n");

  for (i=0; i < lp -> c0; ++i)
    if (lp -> lb [i] != 0. || 
	lp -> ub [i] < .5e30)
      fprintf (f, "%g <= x_%d <= %g\n", 
	       lp -> lb [i], 
	       i,
	       lp -> ub [i]);

  fprintf (f, "\nend\n");

  fclose (f);

  return 0;
}
