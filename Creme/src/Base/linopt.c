/*
 * Name:    linopt.c
 * Author:  Pietro Belotti
 * Purpose: linear optimization (solve one-dimensional MaxFS problems)
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sparse.h"
#include "linopt.h"

#define RTR_BIG_O_CONST 10

#define RTR_THRES_WEAVE 5000 /* up to these lambdas can be sorted
			        without sharing computations */

/*
 * Compare function - used by qsort
 */

inline int compare_abs (const void *a, const void *b) {

  register double x = fabs (* ((const double *) a));
  register double y = fabs (* ((const double *) b));

  if      (x < y) return -1;
  else if (x > y) return  1;
  else            return  0;
}


/*
 * Fast vector operations
 */

#ifdef RTR_USE_PRAGMAS

__inline void calc_lhs_block (double *z, double *coe, double *x, int *pos, register int n) {

#pragma disjoint (*z, *coe, *x, *pos)

  for (; n >= 0; --n)
    *z += coe [n] * x [pos [n]];
}

__inline void calc_lhs (double *z, int n, double *coe, double *x, int *pos) {

  register int j = n / CONCURRENT_FLOPS;

  for (; j > 0; --j, pos += CONCURRENT_FLOPS,
	             coe += CONCURRENT_FLOPS)
    calc_lhs_block (z, coe, x, pos, CONCURRENT_FLOPS);

  calc_lhs_block   (z, coe, x, pos, n % CONCURRENT_FLOPS);
}

#endif


/*
 * find the best point over the segment [x, x + alpha dx]
 */

double scan_segment (double *lambda, int nl, int *dSat, int *fSat) {

  double z = 0;
  int best;

  register double *pl;
  register int i,j;

  for (j  = best = 0,
       i  = nl-1,
       pl = lambda; i >= 0; --i, pl++) {

    if (*pl > 0) j++;
    else         j--;

    if (j>best) {

      best = j;
      if (i) z = 0.5 * (fabs (*pl) + fabs (*(pl+1)));
      else   z = fabs (*pl) + 10 * EPSILON;
    }
    else if (!fSat && (j < -i)) break;
  }

  if (dSat) *dSat = best;
  if (fSat) *fSat = j;

  return z;
}


/*
 * Merge sorted subvectors
 */

double scan_sorted_subvectors (int k, register int nl,
			       double *lambdas, int *displs, int *nlks,
			       int *dSat, int *fSat) {

  if ((nl <= RTR_THRES_WEAVE) ||
      (k >= RTR_BIG_O_CONST * log (nl) / log (2))) {

    /*
     * lambda is short, sort it rather than merge it
     */

    qsort (lambdas, nl, sizeof (double), compare_abs);

    return scan_segment (lambdas, nl, dSat, fSat);
  }
  else {

    int i, j = 0; /* constraints "entered" - constraints "left" */
    int maxj = -1;

    double *best = NULL;
    double minl, minl2, z = 0.;
    int kbest = -1;

    double l2;

    register double l;
    register double **pl = (double **) malloc (k * sizeof (double *));

    for (i=k; i>0; i--) *pl++ = lambdas + *displs++;

    displs -= k;
    pl     -= k;

    /*
     * manual merge
     */

    for (; nl>0; nl--) {

      minl2 = minl = 1e30;

      for (i=0; i<k; i++, pl++, nlks++) if (*nlks) {

	l = fabs (**pl);

	if (*nlks) {

	  if (l < minl) {

	    if ((*nlks>1) && ((l2 = fabs(*(*pl+1))) < minl2)) minl2 = l2;

	    if (minl < minl2) minl2 = minl;

	    minl  = l;
	    best  = *pl;
	    kbest = i;
	  }
	  else if (l < minl2) minl2 = l;
	}
      }

      if (minl2 > 1e9) minl2 = minl + 10 * EPSILON;

      nlks -= k;
      pl   -= k;

      nlks [kbest] --;
      pl   [kbest] ++;

      if (*best < 0) {

	if ((--j < - nl) && !fSat)
	  return z;
      }
      else {

	if (++j > maxj) {

	  maxj = j;

	  if (i) z = 0.5 * (minl + minl2);
	  else   z = minl + 10 * EPSILON;
	}
      }
    }

    free (pl);

    if (dSat) *dSat = maxj;
    if (fSat) *fSat = j;

    return z;
  }
}


/*
 * unidimensional optimization
 */

double one_opt (sparseLP *lp,
		double *dx,
		double *b_Ax,
		double stretch) {

  double *lambda;  /* points of intersection over the segment */

  int nlk = 0, nl, /* number of constraints across the segment (local, global) */
    r, c, i, j, k;

  double **ic, *pc;
  int    **ip, *pp;
  int     *il;
  double  *ub, *lb;

  double z, sum;
  int    *nlks;
  int    *displs;

  r  = lp -> rk;
  c  = lp -> ck;

  ic = lp -> ic;
  ip = lp -> ip;
  il = lp -> il;

  ub = lp -> ub;
  lb = lp -> lb;

  lambda = (double *) malloc ((r+1) * sizeof (double));

  /*
   * build lambda
   */

  for (i=0; i<r; i++) {

    pc = *(ic++);
    pp = *(ip++);

    sum = 0;

#ifdef RTR_USE_PRAGMAS
    calc_lhs (&sum, *il++, pc, dx, pp);
#else
    for (j=*(il++); j>0; j--) sum += (*(pc++) * dx [*(pp++)]);
#endif

    /* does constraint i cross the segment?
     *
     * yes, if Ax < b and A (x + dx) >= b or vice versa
     * and the crossing is before the end of the segment
     */

    z = *b_Ax++;

    if ((((z <= 0) && (sum < 0)) ||
         ((z >  0) && (sum > 0)))
	&& ((z /= sum) < stretch)) {     /* z & sum have the same sign */

      *lambda++ = (sum >= 0) ? z : -z;
      nlk++;
    }
  }

  k = lp->ncpus;

  lambda -= nlk;

  nlks   = (int *) malloc (k * sizeof (int));
  displs = (int *) malloc (k * sizeof (int));

  nl = nlk;

  if (!nl) {
    z=0;
    goto end;
  }

  if (nlk > 1) qsort (lambda, nlk, sizeof (double), compare_abs);
  if (nlk < (int) (r * 0.99)) lambda = (double *) realloc (lambda, (nlk+1) * sizeof (double));

  z = scan_segment (lambda, nlk, NULL, NULL);

 end:

  free (displs);
  free (nlks);
  free (lambda);

  return z;
}
