/*
 * Name:    move.c
 * Author:  Pietro Belotti
 * Purpose: move to next iterate and update all sparse data structures
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "sparse.h"
#include "linopt.h"
#include "rtr.h"

#ifdef RTR_USE_PRAGMAS

/*
 * fast vector operations (use #pragmas)
 */

__inline void build_dxk_block (int *ip, double *ic, double *dxk, double z) {

  register int i = CONCURRENT_FLOPS;

#pragma disjoint (*dxk, *ip, z, *ic)

  for (; i>0; --i) dxk [ip [i]] += z * ic [i];  /* update variation of variable *ip */
}


__inline void build_dxk (int k, int *ip, double *ic, double *dxk, double z) {

  register int i = k / CONCURRENT_FLOPS;

  k %= CONCURRENT_FLOPS;

  for (; i>0; --i) {

    build_dxk_block (ip, ic, dxk, z);

    ip += 10;
    ic += 10;
  }

  /* there is something fishy here... why not call it with k%...? */
}

#endif


/*
 * Update points and aux. vectors for a delta vector
 */

int update (sparseLP *lp, double *x, double *dx, double *b_Ax, 
	    char *sat, double *sum, double maxMod) {

  int i, j, dSatd = 0;

  int     c    = lp->c0;
  int    *vl   = lp->vl;

  int    **colP = lp->vp;
  double **colC = lp->vc;

  double *ub = lp->ub;
  double *lb = lp->lb;

  int    *vp;
  double *vc;

  double z;

  for (i=0; i<c; i++, vl++, colC++, colP++, dx++, x++, ub++, lb++) {

    if ((*dx > EPSILON) || (*dx < EPSILON)) {

      vp = *colP;
      vc = *colC;

      z = *x;

      if      ((*x += maxMod * *dx) > *ub) {*x = *ub;}
      else if  (*x                  < *lb) {*x = *lb;}

      z -= *x;

      /* for all constraints containing this variable, update b_Ax */
      /* and respective indicator                                  */

#ifdef RTR_USE_PRAGMAS
#pragma execution_frequency(very_high)
#endif

      for (j = *vl; j>0 ; --j) {

	register double *pV    = b_Ax + *vp;
	register char   *pS    = sat  + *vp++;
	register double  delta = z    * *vc++;

	if ((*pV) > 0) {
	    
	  if ((delta > 0) || ((*pV) > - delta))
	    *sum += delta;
	  else {
	    *sum -= *pV;
	    if (!*pS) { dSatd++; *pS = SATD; } /* constraint is now fulfilled */
	  }
	}
	else {
	  
	  if (delta > - (*pV)) {
	    *sum += (*pV + delta);
	    if  (*pS) { dSatd--; *pS = UNSATD; } /* constraint is now violated */
	  }
	}

	*pV += delta; /* update b_Ax */
      }

      *dx = 0;
    }
  }

  return dSatd;
}


/*
 * create move direction and move to new x accordingly
 */

int move (sparseLP *lp,     /* LP data  */
	  double   *x,      /* current point  */
	  double   *b_Ax,   /* current violation  */
	  int      *block,  /* array with indices of selected constraints  */
	  char     *sat,    /* indicator function of fulfilled constraints  */
	  double    temp,   /* temperature  */
	  double   *sum     /* total violation (updated here)  */
	  ) {

  static double *dx0 = NULL; /* move direction (global)  */
  static double *dxk = NULL; /* move direction (local)  */

  static char first  = 1;    /* is this the first call to the function?  */

  int i, c, r;

  double *ub, *lb;

  double z, maxMod;

  /* free static malloc'ed arrays  */

  if (!lp) {

    free (dx0);
    free (dxk);
    return 0;
  }

  /* first call, initialize local arrays  */

  c = lp->c0;
  r = lp->rk;

  if (first) { 

    first = 0;

    dx0 = (double *) malloc (c * sizeof (double));
    dxk = (double *) malloc (c * sizeof (double));

    for (i=c; i>0; i--) *dxk++ = 0.0;

    dxk -= c;
  }

  /* Compute dx based on b_Ax  */

  lb   = lp -> lb;
  ub   = lp -> ub;

  for (; *block >= 0; block++) {

    int    *ip = (lp->ip) [*block];
    double *ic = (lp->ic) [*block];

    z  = b_Ax     [*block];

    if (z<0) z = 0;
    else     z = exp (- z / temp);

#ifdef RTR_USE_PRAGMAS
    build_dxk ((lp->il) [*block], ip, ic, dxk, z);
#else
    for (i = (lp->il) [*block]; i>0; i--, ip++)
      dxk [*ip] += z * *ic++; /* update variation of variable *ip   */
#endif
  }

  for (i=c; i>0; --i) *dx0++ = *dxk++;
  dxk -= c;
  dx0 -= c;

  for (i=c; i>0; --i) *dxk++ = 0;
  dxk -= c;

  /* Project back into bounding box  */

  if (lp -> onedim) maxMod = lp -> stretch;
  else              maxMod = lp -> alpha;

  /*
    {
    double mind = 1e40, maxd = -1e40, norm=0;

    for (i=c;i>0;i--, dx0++) {
    if (fabs (*dx0) > maxd) maxd = *dx0;
    if (fabs (*dx0) < mind) mind = *dx0;
    norm += *dx0 * *dx0; 
    }
    norm = sqrt (norm);
    dx0 -= c;
    printf ("||dx0|| = %12.8f (%12.8f,%12.8f); ", norm, mind, maxd);
    }
  */

  for (i=c; i>0; --i, dx0++, x++, ub++, lb++) {

    if      (*dx0 < - EPSILON) {       /* move down, check if below lower bound  */
      if (*x <= *lb) *dx0 = 0;         /* already on boundary   */
      else if ((z = (*lb - *x) / *dx0) < maxMod) maxMod = z;
    }

    else if (*dx0 >   EPSILON) {        /* move up,  check if above upper bound  */
      if (*x >= *ub) *dx0 = 0;          /* already on boundary   */
      else if ((z = (*ub - *x) / *dx0) < maxMod) maxMod = z;
    }

    /*
     * Now maxMod contains a scaling factor for (x + maxMod * dx) to
     * stay within bounding box
     */
  }
  
  dx0 -= c;
  x   -= c;
  ub  -= c;
  lb  -= c;
  
  /* If one-dimensional optimization is used, tune the next point so
   * as to get the most fulfilling point in the segment
   */

  /*  printf ("maxMod = %12.8f ", maxMod);  */


  if (lp -> onedim) { /* create sparse version of dx0  */

    z = one_opt (lp, dx0, b_Ax, maxMod);

    if (z > 0) maxMod = z;
  }

  /*  printf ("-> %12.8f\n", maxMod);  */

  /* {
     double summ = 0;
     for (i=0; i<c; i++, dx0++) summ += *dx0 * *dx0;
     printf ("||dx|| = %.6f\n", log (1e-15+fabs(maxMod * sqrt(summ))) / log (10)); 
     dx0 -= c;
     } */

  /*  Update b_Ax  */

  return update (lp, x, dx0, b_Ax, sat, sum, maxMod); 
}
