/*
 * Name:    rtr.c
 * Author:  Pietro Belotti
 * Purpose: find a maximum feasible subsystem of an infeasible LP
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 */

#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "sparse.h"
#include "init.h"
#include "rtr.h"
#include "move.h"
#include "locsrch.h"
#include "linopt.h"
#include "chooseblock.h"
#include "isfeas.h"

#define MAX_NIMPROV 2000

#define MIN_TEMP    0.0001
#define MAX_TEMP   1e-4

#define INIT_MOMENTUM 40

extern int interrupt;

/*
 * Randomized Thermal Relaxation algorithm
 */

int rtr (sparseLP *lp, /* description of infeasible LP  */
	 char *sol     /* contains feasible subsystem   */
	 ) {

  char lincool = lp -> lincool,
       invcool = lp -> invcool;

  int maxIter  = lp -> nIter;

  double alpha     = lp -> alpha,
         beta      = lp -> beta,
         gammaRate = lp -> gammaRate,
         muRate    = lp -> muRate;

  int      nSatd = 0,
       nSatd_loc = 0,
      best_nSatd = 0,
      nImprov = 0,
      nIter;

  int dnLoc, /* improvement with variable local search  */
      dnAms; /* improvement with rtr  */

  int rtrMomentum = INIT_MOMENTUM;

  char whichmove = USE_RTR;

  double *x;        /* current solution  */
  double *b_Ax;     /* current ineq. violation (>0 if violating)  */
  char   *satd;     /* satd [i] = 1 if i-th ineq. is satisfied, 0 otherwise  */
  int    *block;    /* set of inequalities used for computing x (k+1)  */

  double temperature = 0, gamma = 1, sumViol, sumViol_loc, mu = 1;

  double time = 0;

  if (lp -> my_id == 0) time = CoinCpuTime();

  if (!(x    =(double*) malloc ((lp->c0)   * sizeof(double)))) fprintf (stderr,"Error: malloc failed.\n");
  if (!(b_Ax =(double*) malloc ((lp->rk)   * sizeof(double)))) fprintf (stderr,"Error: malloc failed.\n");
  if (!(satd =(char  *) malloc ((lp->rk)   * sizeof(char))))   fprintf (stderr,"Error: malloc failed.\n");
  if (!(block=(int   *) malloc ((lp->rk+1) * sizeof(int))))    fprintf (stderr,"Error: malloc failed.\n");

  /*
   *   Random initialization of all variables
   */

  if (lp -> my_id == 0) {

    printf ("   #iter     |mfs|        Temp.      time   |bl.|   %%covg.\n");
    printf ("==========================================================\n");
  }

  for (nIter = 0; 

       ((lp->timelimit < 0) || (CoinCpuTime() - time <= lp -> timelimit)) &&
	 (nIter < maxIter) && 
	 !interrupt && 
	 (nSatd < lp->r0); 

       nIter++) {

    /*
     * (Re)start from randomized solution
     */

    if ((!(nIter % (lp->restFreq)) ||    /* time to restart  */
 	 (nImprov > MAX_NIMPROV)   ||    /* too many non-improving iterations  */
	 (temperature*gamma<MIN_TEMP))   /* temperature too low  */
	&& (nSatd <= best_nSatd)) {      /* aspiration criteria: do not restart   */
                                         /* if good solution found  */

      nImprov = 0;

      /*      if ((anneal *= 2.0) > 0.4) anneal = 0.4;  */

      init_x (lp, x);

      /*
       *  satd & b_Ax initialized using current x values.
       */

      nSatd_loc = init_sat (lp, satd, b_Ax, x, &sumViol_loc);

      sumViol = sumViol_loc;

      /*
       * TO DO: Should alpha be computed locally?
       *        Should we use sumViol_loc instead?
       */ 

      /*
       *  CAUTION! mu decreases in non improving iterations, but if
       *  other rules are used the alpha needs be recomputed globally
       */

      temperature = alpha * sumViol / (lp -> rk - nSatd_loc); 

      if (temperature > MAX_TEMP) 
	temperature = MAX_TEMP;

      gamma = 1;
      mu    = 0.3;
    }

/*    if ((maxIter<1000) || !(nIter % (maxIter/1000)))  */
/*      fprintf (stderr, "\r%.1f%%", (nIter+1) * 100.0 / maxIter);  */

    nSatd = nSatd_loc;

    /*
     *   Found a better point?
     */

    if (nSatd > best_nSatd) {

      int i, feas = 1;

      int feask = (isFeas (lp, satd, x, &nSatd_loc)) == 3;

      feas = feask;

      if (lp -> my_id == 0) {

	printf ("%8d%c %9d %11.2f %8.2f %9d %7.3f %c\n",
		nIter, ((whichmove == USE_LOCSRCH)? '+' : ' '), 
		nSatd, temperature * gamma,
		CoinCpuTime () - time,
		mymax (1, (int) (mu * (lp->r0 - nSatd))),
		(int) (1000.0 * (100.0 *  nSatd / lp -> r0)) / 1000.0,
		feas? '*' : '!');

	fflush (stdout);

	/*
	 * test actual feasibility of the subsystem found
	 */

	/*	fprintf (stderr, "\n\r%.1f%%", (nIter+1)*100.0/maxIter); fflush (stderr);  */
      }

      best_nSatd = nSatd;

      /*
       * save current solution
       */

      /* #pragma disjoint  */

      for (i=lp->rk-1; i>=0; --i) *sol++ = *satd++;
      satd -= lp->rk;
      sol  -= lp->rk;

      nImprov = 0;

      /*
       * this is a good point, increase block percentage
       */

      if ((mu *= 1.1) > 1) mu = 1;

      rtrMomentum  = INIT_MOMENTUM;
    }

    else {
                                           /* no improvement:  */
      if (rtrMomentum > 0) --rtrMomentum;  /* - decrease importance of rtr vs. variable local search  */
      mu /= muRate;                        /* - shrink block fraction  */
      nImprov++;                           /* - increase no. non-improving iterations  */
    }

    /*
     * Main loop
     */

    if (whichmove == USE_RTR) { /* primal move (pure rtr)  */

      /*
       * select a set of inequalities to be used in computing x(k+1).
       */

      choose_block (lp, block, satd, nSatd_loc, 
		    mymax (1, (int) (mu * (lp -> rk - nSatd_loc))),
		    b_Ax, sumViol);

      /*
       * calc. x(k+1), update b_Ax, satd
       */

      nSatd_loc += (dnAms = move (lp, x, b_Ax, block, satd, temperature * gamma, &sumViol));

      if ((dnAms <= 0) && (lp -> locsea) && (rtrMomentum == 0)) {

	/*
	 * too many non-improving iterations with rtr, next iteration
	 * try variable local search
	 */

        whichmove    = USE_LOCSRCH;
        rtrMomentum  = INIT_MOMENTUM;
      }
    }

    else {      /* dual move: variable local search  */

      nSatd += (dnLoc = locsrch (lp, x, b_Ax, satd, &sumViol));

      if (dnLoc <= 0) whichmove = USE_RTR; /* no improvement with variable local search  */
                                           /* choose rtr as next method  */
    }

    /*
     * End of main loop
     */

    if      (lincool) gamma -= 1.0 / maxIter;         /* update gamma linearly (temperature scaling factor)  */
    else if (invcool) gamma *= (nIter+1) / (nIter+2); /* update gamma inverse-linearly  */
    else              gamma /= gammaRate;             /*              geometrically  */

    /* update temperature:   */
    /* convex combination (beta) of previous temperature and  */
    /* current total violation averaged on violated constraints  */

    if (nSatd_loc < lp -> rk)
      temperature =      beta  * temperature +
                    (1 - beta) * alpha * sumViol / (lp -> rk - nSatd_loc);
  }

  if (lp -> my_id == 0) {

    if      (interrupt)       printf ("User interrupt\n");
    else if (nSatd >= lp->r0) printf ("All ineqs satisfied\n");
    else                      printf ("Completed\n");
  }

  if (lp -> my_id == 0) printf ("Total time: %.2f\n", CoinCpuTime() - time);

  /*
   * check actual feasibility
   */

  switch (isFeas (lp, satd, x, &nSatd_loc)) { /* 3 indicates the mfs is actually feasible  */

    case -1: fprintf (stderr, "\rBounds Violated\n");                               break;
    case  0: fprintf (stderr, "\rError: mfs with no sense\n");                      break;
    case  1: fprintf (stderr, "\rWarning: mfs includes unsatisfied ineqs\n");       break;
    case  2: fprintf (stderr, "\rWarning: mfs does not include satisfied ineqs\n"); break;
    case  3: if (nSatd == lp->r0) fprintf (stderr, "\rProblem Feasible\n");         break;
    default: fprintf (stderr, "\rCase out of bounds\n");
  }

  /*
   * free static malloc'ed pointers in subroutines
   */

  move    (NULL, NULL, NULL, NULL, NULL, 0, NULL);
  locsrch (NULL, NULL, NULL, NULL, NULL);

  free (b_Ax);
  free (satd);
  free (block);
  free (x);

  return best_nSatd;
}
