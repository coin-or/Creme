/*
 * Name:    sparse.h
 * Author:  Pietro Belotti
 * Purpose: data structures for sparse infeasible LP
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 */

#ifndef SPARSE_H
#define SPARSE_H

#define COEFF_TYPE double

#define SATD   1
#define UNSATD 0

#define WEIGHT_STEP 100
#define WEIGHT_MAX  10000

#define EPSILON 1e-8

#define CONCURRENT_FLOPS 10

#define PRINTF(x) {printf ("%2i: ", lp->my_id); printf (x);}
#define mymax(a,b) ((a) < (b) ? (b) : (a))
#define mymin(a,b) ((a) > (b) ? (b) : (a))

/*
 *    Sparse LP
 *
 *    contains a sparse description of the coefficient matrix
 *    (specified by pairs (index, value) and parameters for the rtr
 *    algorithm
 *
 */

typedef struct _sparseLP {

  int       r0,c0; /* no. rows/columns  */
  int       rk,ck; /* no. rows/columns in this chunk  */
  int       nnzk;  /* # nonzero  */

  int         *il; /* inequality length	(nonzero coefficient of j-th row)  */
  COEFF_TYPE **ic; /* inequality coefficient a_j  */
  int        **ip; /* inequality index of a_j (i.e., j)  */

  COEFF_TYPE *rhs; /* right hand side b of ax >= b  */
  COEFF_TYPE *rlb; /* constraint's lower bound (the b in "b <= ax <= c"), will replace rhs  */
  COEFF_TYPE *rub; /* constraint's upper bound (the c in "b <= ax <= c"), will replace rhs  */

  int         *vl; /* variable length (nonzero coefficient of i-th column)  */
  COEFF_TYPE **vc; /* variable coefficient a_j  */
  int        **vp; /* variable index of a_j (i.e., j)  */

  COEFF_TYPE  *lb; /* mandatory lower bound for variables  */
  COEFF_TYPE  *ub; /*           upper  */

  int ntaut;       /* number of tautologies (constraints fulfilled by any x in [l,u])  */
  int niis;        /* number of 0-iis       (            violated                   )  */

  /*  double *cum_weights;  // cumulated weights  */
  char *chosen;         /* 1 if constraint is included in the block, 0 otherwise  */

  char  noprep,    /* no preprocessing  */
    bigm,          /* only write big-M file  */
    onedim,        /* one-dimensional search  */
    locsea,        /* variable local search  */
    lincool,       /* linear temperature decrease  */
    norm,          /* normalize solution  */
    dblrand,       /* use double randomization  */
    invcool;       /* inverse linear temperature decrease  */

  int   nIter,     /* number of iterations  */
    restFreq,      /* restart every this iterations  */
    blkcard,       /* block cardinality for variable local search (see locsea)  */
    ncpus,         /* number of parallel threads  */
    my_id;         /* process identifier  */

  double alpha,    /* scaling factor of average violation  */
    beta,          /* convex combination parameter on segment (old temperature, new temperature)  */
    gammaRate,     /* scaling value of the temperature (decreases either linearly or exponential)  */
    muRate,        /* controls decrease of block size after non-improving iterations  */
    timelimit,     /* maximum cpu time  */
    infinity,      /* infinity  */
    stretch;       /* elongation of dx0 in one-dimensional optimization  */

} sparseLP;

void userInterrupt ();

#endif
