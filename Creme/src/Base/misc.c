/*
 * Name:    misc.c
 * Author:  Pietro Belotti
 * Purpose: miscellaneous functions for reading BZ files, printing LPs, etc
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


/*
 * Print LP (for debugging purposes)
 */

void printLP (int *mfs, sparseLP *lp) {

  int i,j,k,l;

  printf ("Minimize\n obj: ");

  k=0;
  for (i=0; i<lp->c0; i++) {
    if (!(++k % 9))
      printf ("\n    ");
    printf (" + 1 x%d", i);
  }

  printf("\nSubject To\n");

  for (i=0; (!mfs && (i<lp->r0)) || (mfs && (mfs[i] >= 0)); i++) {

    l = mfs ? (mfs[i]) : i;

    printf ("c%d:", l);
    k=0;

    for (j=0; j < lp -> il [l]; j++) {
      if (!(++k % 5)) {printf ("\n    "); k++;}
      printf ("%c%f x%d ", (((lp -> ic) [l] [j]) > 0 ? '+' : ' '),
                             (lp -> ic) [l] [j],
                             (lp -> ip) [l] [j]);
    }

    if (!(k++ % 5)) printf ("\n");
    printf (">= %f\n", lp -> rhs [l]);
  }

  printf ("Bounds\n");

  for (i=0; i<lp->c0; i++)
    printf ("%f <= x%d <= %f\n", lp->lb [i], i, lp->ub [i]);

  printf ("End\n");
}


/*
 * Remove LP from memory
 */

void clearLP (sparseLP *lp) {

  int i;

  for (i=lp->c0 - 1; i>=0; i--) free (lp -> vc [i]);
  for (i=lp->c0 - 1; i>=0; i--) free (lp -> vp [i]);

  free (lp->vl);

  free (lp->vp);
  free (lp->vc);

  free (lp->ic);
  free (lp->ip);
  free (lp->il);

  free (lp->rhs);

  free (lp->chosen);
}

/****************************************************/
/*                                                  */
/*  Functions to read from a bzip compressed file   */
/*                                                  */
/****************************************************/

/*
 *  get char from .bz2 file
 *
 * ONLY FOR DEBUGGING PURPOSES
 * DO NOT USE WITH bzGetDouble or bzGetInt
 */

int bzgetchar (BZFILE *f, char *c) {

  static char str [MAX_STR];
  static char *p = NULL;
  static int count = 0;

  int status, i;

  if (!p) {

    i = BZ2_bzRead (&status, f, str, MAX_STR);
    p = str;

    if (i < MAX_STR) str [i] = 1;
    count = 0;
  }

  *c = *p++;

  if (*c == 1) return 1;

  if (++count == MAX_STR) p = NULL;

  return 0;
}


/*
 *  get double from .bz2 file
 */

int bzgetdbl (BZFILE *f, double *val) {

  static char str [2 * MAX_STR];
  static char *p = NULL;

  int i, status;

  if (!p) {

    i = BZ2_bzRead (&status, f, str, 2 * MAX_STR);
    p = str;

    if (i < 2 * MAX_STR) str [i] = 1;
  }

#ifdef RTR_USE_PRAGMAS
#pragma execution_frequency(very_high)
#endif
  {

    register int i,j;

    char flag = 0;

    i = strcspn (p,   "+-0123456789eE.");

    if (p [i] == 1) {
      p = NULL;
      return 1;
    }

    j = strspn  (p+i, "+-0123456789eE.");

    if (p [i+j] == 1) flag = 1;

    p [i+j] = 0;

    *val = atof (p+i);

    if (flag || (*(p += i+j+1) == 1)) {
      p = NULL;
      return 1;
    }

    if (p >= str + MAX_STR) {

      p -= MAX_STR;
      memcpy (str, str + MAX_STR, MAX_STR);
      i = BZ2_bzRead (&status, f, str + MAX_STR, MAX_STR);
      if (i<MAX_STR) str [MAX_STR + i] = 1;
    }
  }

  return 0;
}


///////////////////////////
//
// get integer from .bz2 file
//

inline int bzgetint (BZFILE *b, register int *val) {

  double x;
  register int ret_val;

  ret_val = bzgetdbl (b, &x);

  *val = floor (x);

  return ret_val;
}


/////////////////////////////////////////////////////////
//
// Auxiliary functions to reallocate vectors of double/int
//

inline void reallocate_double (int req, int *size, double **buf) {

  if (req >= *size)
    *buf = (double *) realloc (*buf, (*size = (MALLOC_BLOCK * (req/MALLOC_BLOCK + 1))) * sizeof (double));
}

inline void reallocate_int (int req, int *size, int **buf) {

  if (req >= *size)
    *buf = (int *) realloc (*buf, (*size = (MALLOC_BLOCK * (req/MALLOC_BLOCK + 1))) * sizeof (int));
}


/*
 * Transpose coefficient matrix
 */

void create_transpose (sparseLP *lp) {

  int i, k;

  int c = lp->c0;
  int r = lp->rk;

  int     *vl = lp -> vl = (int     *) malloc (c * sizeof (int));
  int    **vp = lp -> vp = (int    **) malloc (c * sizeof (int    *));
  double **vc = lp -> vc = (double **) malloc (c * sizeof (double *));

  int     *il = lp -> il;
  int    **ip = lp -> ip;
  double **ic = lp -> ic;

  for (i=c; i>0; i--) {

    *vl++ = 0;
    *vc++ = NULL;
    *vp++ = NULL;
  }

  vl -= c;
  vc -= c;
  vp -= c;

  for (i=0; i<r; i++) {

    register int    *pp = *ip++;
    register double *pc = *ic++;

    for (k=*il++; k>0; k--, pc++, pp++) {

      register int j = vl [*pp] ++;

      if (!(j % MALLOC_BLOCK)) {

	vc [*pp] = (double *) realloc (vc [*pp], (j + MALLOC_BLOCK) * sizeof (double));
	vp [*pp] = (int    *) realloc (vp [*pp], (j + MALLOC_BLOCK) * sizeof (int));

	if (!(vp[*pp]) || !(vc[*pp])) {
	  perror ("create_transpose");
	  exit (-1);
	}
      }

      vc [*pp] [j] = *pc;
      vp [*pp] [j] = i;
    }
  }


  for (i=c; i>0; i--, vl++, vp++, vc++) {

    *vp = (int    *) realloc (*vp, *vl * sizeof (int));
    *vc = (double *) realloc (*vc, *vl * sizeof (double));
  }
}


//////////////////////////////////////
//
// Get the norm of a vector of size l
//

static inline double get_norm (register double *v, register int l) {

  register double norm = 0;
  for (; l>0; --l, ++v) norm += *v * *v;
  return sqrt (norm);
}
