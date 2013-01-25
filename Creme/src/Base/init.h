/*
 * Name:    init.h
 * Author:  Pietro Belotti
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#ifndef INIT_H
#define INIT_H

#include "sparse.h"

void init_x   (sparseLP *lp, double *x);

int  init_sat (sparseLP *lp, char   *sat, double *b_Ax, double *x, double *temp);


/*
 *   Initialize b-Ax [j], sat [j], scalars sum_of_violation return
 *   number of constraints fulfilled
 *
 *   It is run for each constraint, sat and b_Ax are the current
 *   positions in the respective arrays
 *
 *   Parallel version: run by all slaves on own chunk; send returned
 *   value to master
 */

#ifdef RTR_USE_PRAGMAS

/*
 * single block operation
 */

static __inline void calc_lhs_block (double *z, double *coe, double *x, int *pos, register int n) {

#pragma disjoint (*z, *coe, *x, *pos)

  for (; n >= 0; --n) 
    *z += coe [n] * x [pos [n]];
}

/*
 * divide in blocks of CONCURRENT_FLOPS size and use pragmas in calc_lhs_block
 */

static __inline void calc_lhs (double *z, int n, double *coe, double *x, int *pos) {

  register int j = n / CONCURRENT_FLOPS;

  for (; j > 0; --j, pos += CONCURRENT_FLOPS, 
	             coe += CONCURRENT_FLOPS)
    calc_lhs_block (z, coe, x, pos, CONCURRENT_FLOPS);

  calc_lhs_block   (z, coe, x, pos, n % CONCURRENT_FLOPS);
}
#endif

#endif
