/*
 * Name:    init.c
 * Author:  Pietro Belotti
 * Purpose: initialize data structures of the sparse LP
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "sparse.h"
#include "init.h"

/*
 * Initialize variables
 */

void init_x (sparseLP *lp, double *x) {

  double *lb = lp -> lb,
         *ub = lp -> ub;

  if (lp -> my_id == 0) {

    register int i;

    for (i = lp->c0; i > 0; --i, ++ub, ++lb)

      if (*ub - *lb > 1e6)
	if (fabs (*lb) < 1e5) *x++ = *lb + drand48 ();
	else                  *x++ = 0;/* *ub - drand48 (); */
      else                    *x++ = *lb + drand48 () * (*ub - *lb);

    x  -= (lp -> c0);
    lb -= (lp -> c0);
    ub -= (lp -> c0);
  }
}


/*
 * Initialize auxiliary vectors based on x
 */

int init_sat (sparseLP *lp,   /* LP data                                                     */
	      char *sat,      /* sat [j] is 1 if j-th constraint fulfilled, 0 otherwise      */
	      double *b_Ax,   /* b_Ax [j] is the (positive) violation of the j-th constraint */
	      double *x,      /* initial solution                                            */
	      double *sumViol /* sum of all violations                                       */
	      ) {

  register int i;

  int      ns = 0;

  double **ic  = lp -> ic, *coe;
  int    **ip  = lp -> ip, *pos;

  double  *rhs = lp -> rhs;
  int     *il  = lp -> il;

  double z;

  *sumViol = 0;

  for (i=lp->rk; i--;) {

    coe =   *  ic ++;
    pos =   *  ip ++;
    z   = - * rhs ++;

#ifdef RTR_USE_PRAGMAS
    calc_lhs (&z, *il++, coe, x, pos);
#else
    {
      register int j;
      for (j = *il++; j>0; j--)
	z += *coe++ * x [*pos++];
    }
#endif

    *b_Ax++ = -z;

    if (z<0) {

      *sat++ = UNSATD;
      *sumViol -= z;
    }
    else {
      *sat++ = SATD;
      ns ++;
    }
  }

  return ns;
}
