/*
 * Name:    locsrch.c
 * Author:  Pietro Belotti
 * Purpose: search for best point by fixing all variables except one
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include "sparse.h"
#include "linopt.h"
#include "locsrch.h"

#define EPS 1e-6
#define ALPHA 0.1

int compare_frontier (const void *one, const void *two) {

  if (((const frontier *) one) -> thres < ((const frontier *) two) -> thres) return -1;
  if (((const frontier *) one) -> thres > ((const frontier *) two) -> thres) return  1;
  return 0;
}

int compare_dvar (const void *one, const void *two) {

  if (((const dvar *) one) -> nsi > ((const dvar *) two) -> nsi) return -1;
  if (((const dvar *) one) -> nsi < ((const dvar *) two) -> nsi) return  1;
  return 0;
}

int locsrch2 (sparseLP *lp,
	      /*double   *x, */
	      double   *b_Ax
	      /*double   *dx */
	      /*	     char     *satd, */
	      /*     double   *sumViol */
	     ) {

  /*
  *  for each var {
  *
  *    1) get lambdas, pluses, minuses;
  *
  *    2a) sort them according to [alpha * pluses + (1-alpha) *
  *    pluses/(pluses+minuses)] select J:|J|=lp->blkcard, first |J|
  *    in ordering
  *
  *    or
  *
  *    2b) select randomly with probability dependent on [alpha *
  *    pluses + (1-alpha) * pluses/(pluses+minuses)] or uniformly
  *
  *    3) call one_opt on convex combination
  *
  *  }
  */


  int r;
  int c;

  int i, j;

  static double *impr = NULL;
  static double *avg  = NULL;

  double **vc = NULL;
  int    **vp = NULL;
  int     *vl = NULL;

  if (!lp) {

    if (impr) free (impr);
    if (avg)  free (avg);

    return 0;
  }

  r = lp -> rk;
  c = lp -> c0;

  if (!impr) {

    impr = (double *) malloc (c * sizeof (double));
    avg  = (double *) malloc (c * sizeof (double));
  }

  /*
   * Roughly check variables that can have good improvement
   */

  for (i=0; i<c; i++) {

    int up=0, down=0, bad=0;

    register int    *pp = *vp++;
    register double *pc = *vc++;

    int l = *vl++;
    int maximp;

    for (j=l; j>0; j--, pp++, pc++) {

      register double viol = b_Ax [*pp];

      if        (viol < - EPSILON) {           /* constraint fulfilled */
	if       (*pc >   EPSILON) bad++;      /* would violate it     */
	else {if (*pc < - EPSILON) down++;}}
      else if   (viol >   EPSILON) {
	if       (*pc < - EPSILON) bad++;
	else {if (*pc >   EPSILON) up++;}}
    }

    maximp = mymax (up, down);

    /*
     * proportion is important, but the total is more
     */

    *impr++ = ALPHA * (maximp / l) + (1-ALPHA) * (maximp / r);

    /*

intersect at

    lb < x+d < ub

    e_j > 0
    b - A (x + d e_j) = b - Ax - d A e_j
b-Ax >0
b-Ax - d a_j <0   -> a_j > 0

b-Ax <0
b-Ax - d a_j >0   -> a_j < 0

up++

 lb < x + d e_j < ub

 (lb - x) < d < (ub - x)

    e_j < 0
    b - A (x + d e_j) = b - Ax - d A e_j

    */

  }

  impr -= c;

  vc -= c;
  vp -= c;

  /*
  for (i=0; i<c; i++) {
    if (drand48 () < *impr++) *dx++ =
  }
  */

  return 0;
}

int locsrch (sparseLP *lp,
	     double   *x,
	     double   *b_Ax,
	     char     *satd,
	     double   *sumViol
	     ) {

  static frontier *cuts   = NULL;  /* threshold  */
  static dvar     *deltas = NULL;  /* delta vars */
  static double   *dx     = NULL;  /* delta x    */

  dvar *deltas_save;

  int i,j,k, max, argmax, ncuts, cursat, ndelta;

  double th, Delta;

  double *coe;
  int    *pos;

  double **vc;
  int    **vp;

  int    *vl;
  double *lb;
  double *ub;

  int r, c;

  double *rhs;

  double *pV;
  char   *pS;

  double step;
  int dSatd;

  if (!lp) {
    if (cuts) {
      free (deltas);
      free (cuts);
      free (dx);
    }
    return 0;
  }

  if (!cuts) {

    /*
     * non-skeleton implementation, rk or r0?
     */

    deltas = (dvar     *) malloc ( (lp -> c0)    * sizeof (dvar));
    cuts   = (frontier *) malloc (((lp -> rk)+1) * sizeof (frontier));
    dx     = (double   *) malloc ( (lp -> c0)    * sizeof (double));
  }

  vc = lp -> vc;
  vp = lp -> vp;

  vl = lp -> vl;
  lb = lp -> lb;
  ub = lp -> ub;

  r  = lp -> rk;
  c  = lp -> c0;

  rhs = lp -> rhs;

  for (i=ndelta=0; i<c; i++, x++, vl++, vp++, vc++, ub++, lb++) if (*vl) {

    coe = *vc;
    pos = *vp;

    /*    printf ("x%d: ", i);*/

    for (ncuts=k=cursat=0, j=*vl; j>0; j--, pos++, coe++) {

      if (((th = b_Ax [*pos] / *coe + *x) > *lb + EPS) && (th < *ub - EPS)) {

	cuts -> thres = th;

        if (*coe > EPS) {cuts -> sense = 'G';      if (*x > th + EPS) cursat ++;}
        else            {cuts -> sense = 'L'; k++; if (*x < th - EPS) cursat ++;}

	/*        printf ("(%d) %c=%.4f ", *pos, (cuts -> sense=='G')?'>':'<', th);*/
	cuts++; ncuts++;
      }
    }

    /*    printf ("\n");*/

    cuts -> thres = *ub;
    cuts -> sense = 'L';       /* cut */

    cuts -= ncuts;

    qsort ((void *) cuts, (size_t) ncuts, (size_t) sizeof (frontier), compare_frontier);

    max    = k;
    argmax = 0;

    for (j=0; j<ncuts; j++, cuts++)

      if      (cuts -> sense == 'L') k--;
      else if (++k>max)              {max = k; argmax = j+1;}

    cuts -= ncuts;

    if (cursat < max) {

      deltas -> nsi = max - cursat;

      if (argmax==0) deltas -> delta = (*lb + cuts -> thres) / 2 - *x;
      else           deltas -> delta = (cuts [argmax].thres + cuts [argmax-1].thres) / 2 - *x;

      deltas -> index = i;
      deltas++; ndelta++;
    }
  }

  deltas -= ndelta;

  vl -= c;
  lb -= c;
  ub -= c;

  x  -= c;

  if (!ndelta) return 0;

  qsort ((void *) deltas, (size_t) ndelta, (size_t) sizeof (dvar), compare_dvar);

  deltas_save = deltas;

  for (j=0; (j<ndelta) && (deltas -> nsi); j++, deltas++) {
    dx [deltas -> index] = deltas -> delta;
  }

  deltas = deltas_save;

  /* Update satd, x, b_Ax */

  if (lp -> blkcard == 1) {

    step = deltas -> delta;
    k    = deltas -> index;

    /*
    printf ("x%d: %d (%.3f%+.3f)\n", deltas [0] . index,
                                     deltas [0] . nsi, x [k],
                                     deltas [0] . delta); fflush (stdout);
    */

    dSatd = 0;

    pos = lp -> vp [k];
    coe = lp -> vc [k];

    if      ((x [k] += step) > ub [k]) {x [k] = ub [k];}
    else if ( x [k]          < lb [k]) {x [k] = lb [k];}

    for (j = vl [k]; j>0 ; j--) {

      pV    = b_Ax  + *pos;
      pS    = satd  + *pos++;
      Delta = - step * *coe++;

      if ((*pV) > 0) {

        if ((Delta > 0) || ((*pV) > - Delta))
          *sumViol += Delta;
        else {
          *sumViol -= *pV;
          if (!*pS) { dSatd++; *pS = SATD; }
	  /*          printf ("con %d now satd", *(pos-1));*/
        }
      }
      else
        if (Delta > - (*pV)) {
	  /*          printf ("con %d now violated\n", *(pos-1));*/
          *sumViol += (*pV + Delta);
          if (*pS) { dSatd--; *pS = UNSATD; }
        }
      *pV += Delta;
    }

    return dSatd;
  }
  else {

    for (j=0; j<c; j++) *dx++ = 0;

    dx -= c;

    deltas_save = deltas;

    for (j=0; (j<ndelta) && (j<lp->blkcard) && (deltas -> nsi); j++, deltas++) {
      /*      printf ("[%5d]: %4d (%.2f)\n", deltas [j] . index, deltas [j] . nsi, deltas [j] . delta);
	      fflush (stdout);*/
      dx [deltas -> index] = deltas -> nsi * deltas -> delta;
    }

    deltas = deltas_save;

    /*
    for (i=0; i<c;) {
      printf ("%3d %.2f ", i, dx [i]);
      if (!(++i % 20)) printf ("\n");
      fflush (stdout);
    }
    printf ("\n");
    */

    if (0) { /* (lp->algorithm == RTR) { */
/*
	dSatd = 0;
	k = 0;

	maxMod = 100.0;

	for (i=0; i<c; i++, dx++, x++, ub++, lb++) if (*tv++) {

	  if      (*dx < - EPSILON) {
	    if (*x < *lb + EPSILON) *dx = 0;
	    else if ((z = (*lb - *x) / *dx) < maxMod) maxMod = z;
	  }
	  else if (*dx >   EPSILON) {
	    if (*x > *ub - EPSILON) *dx = 0;
	    else if ((z = (*ub - *x) / *dx) < maxMod) maxMod = z;
	  }
	  else *dx = 0;

//        if (((*dx>0) && ((z = (*ub - *x) / *dx) < maxMod))  ||
//	      ((*dx<0) && ((z = (*lb - *x) / *dx) < maxMod)))   maxMod = z;

	}

	tv -= c;

	dx -= c;
	x  -= c;

	ub -= c;
	lb -= c;

	if (maxMod <= 0.0) maxMod = 1.0;

	for (i=0; i<c; i++, vl++, colC++, colP++, dx++, x++, ub++, lb++, tv++) {
	  if (*tv)
	  {
	    *tv = 0;

	    vp = *colP;
	    vc = *colC;

	    z = *x;

	    if      ((*x += maxMod * *dx) > *ub) {*x = *ub;}
	    else if (*x          < *lb) {*x = *lb;}

	    z -= *x;

	    for (j = *vl; j>0 ; j--) {

	      pV    = b_Ax + *vp;
	      pS    = sat  + *vp++;
	      delta = z * *vc++;

	      if ((*pV) > 0)
		if ((delta > 0) || ((*pV) > - delta))
		  *sum += delta;
		else {
		  *sum -= *pV;
		  if (!*pS) { dSatd++; *pS = SATD; }
		}
	      else
		if (delta > - (*pV)) {
		  *sum += (*pV + delta);
		  if (*pS) { dSatd--; *pS = UNSATD; }
		}
	      *pV += delta;
	    }
	    *dx = 0;
	  }
	}
	dx -= c;
	tv -= c;

	ub -= c;
	lb -= c;
	x  -= c;
*/
    }
    /*    else return linOpt (lp, x, dx, b_Ax, satd, sumViol, 1/temp); */
  }

  return 0;
}
