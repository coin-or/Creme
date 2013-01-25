/*
 * Name:    isfeas.c
 * Author:  Pietro Belotti
 * Purpose: check actual feasibility of the feasible subsystem found
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdio.h>

#include "sparse.h"

/*
 * feasibility test
 */

int isFeas (sparseLP *lp,   /* sparse LP data                          */
	    char     *sat,  /* alleged fulfilled constraints           */
	    double   *x,    /* current point                           */
            int      *n     /* alleged number of constraints satisfied */
	    )  {

  double sum;
  int fok = 1;
  int uok = 1;
  int nsatd = 0;
  register int i;

  double *rhs = lp->rhs;
  double *lb  = lp->lb,
         *ub  = lp->ub;

  double **rowC = lp->ic;
  double  *ic;

  int    **rowP = lp->ip;
  int     *ip, *il = lp->il;

  for (i=lp->rk; i>0; i--, sat++) {

    sum = - *rhs++;

    ic = *rowC++;
    ip = *rowP++;

#ifdef RTR_USE_PRAGMAS
    calc_lhs (&sum, *il++, ic, x, ip);
#else
    {
      register int j;

      for (j=*il++; j>0; j--)
	sum += *ic++ * x [*ip++];
    }
#endif

    if (sum >= 0) {
      ++nsatd;
      if (!*sat) uok = 0;
    }
    else if (*sat) fok = 0;
  }

  for (i=lp->c0; i>0; i--, lb++, ub++, x++)
    if ((*x < *lb) || (*x > *ub)) {
      /*      printf ("%d: x_%d = %.5f out of bounds [%.4f,%.4f]\n", 
      	      lp->my_id, lp->c0 - i, *x, *lb, *ub);*/
      return -1;
    }

  if (nsatd != *n) {
    printf ("Thread %d: Warning, numbers of satisfied inequalities differ after check (%d != %d)\n", lp->my_id, *n, nsatd);  
    *n = nsatd;
  }

  /*  if (2*fok+uok!=3) printf ("%d: return %d\n", lp->my_id, 2*fok+uok); */

  return 2 * fok + uok;
}
