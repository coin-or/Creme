/*
 * Name:    chooseblock.c
 * Author:  Pietro Belotti
 * Purpose: select subset of violated inequalities for the update direction
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
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

/*
 * block choice routine
 */

void choose_block (sparseLP *lp     /* LP data                                     */
		   ,int *block      /* arrays with row indices to use (last is -1) */
		   ,char *sat       /* fulfilled constraints			   */
		   ,int nSat        /* number of fulfilled constraints 		   */
                   ,int cardinality /* cardinality of the block			   */
		   ,double *b_Ax    /* violation vector				   */
		   ,double sumviol  /* sum of all violations			   */
		   ) {

  register int j=0, k;

  const int r  = lp->rk;
  char* chosen = lp->chosen;

  /*
    double cum_weight = 0.0;
    double* cum_weights = lp->cum_weights;
  */

  if (r == nSat) return;

  if (!(lp->dblrand))

    while (j<cardinality) {

      if (!(sat [k = floor (drand48() * r)]) &&
	  !(chosen [k])) {

	chosen [k]   = 1;
	block  [j++] = k;
      }
    }
  else {

    double relaxed_divider = sumviol;
    int cnt = 0;

    while (j<cardinality) {

      if (!(sat [k = floor (drand48() * r)]) &&
	  !(chosen [k]) &&
	  (drand48 () < exp (- b_Ax [k] / relaxed_divider))) {

	chosen [k]   = 1;
	block  [j++] = k;
      }

      if (cnt++ > r) {

	relaxed_divider *= 10;
	cnt = 0;
      }
    }
  }

  /*

  for (i = 0, j = 0; i < r; ++i) {
  if (!sat[i]) {
  cum_weight += exp(-b_Ax[i]/sumviol);
  }
  cum_weights[i] = cum_weight;
  }

  for (i = 5 * cardinality - 1, j = 0; j < cardinality && i >= 0; --i) {

  const double key = drand48() * cum_weight;
  const int select = dsearch (key, cum_weights, r);
  if (chosen[select] == 0) {

  block[j] = select;
  chosen[select] = 1;
  ++j;
  }
  }
  */

  block [j--] = -1;

  for (; j >= 0; --j) chosen [block [j]] = 0;
}


/*

inline int cum_weight_comp2(const void* w0, const void* w1)
{
	const double diff = (*(double*)w0) - (*(double*)w1);
 	return diff < 0 ? -1 : ( diff > 0 ? 1 : 0);
}

*/

/*
 * TO DO: make it faster!
 */

/*
static inline int lookup (const double key,
			  const double *v,
			  const int length) {

  register int i;

  for (i=0; i<length; i++) if (v [i] >= key) return i;

  return length;
}
*/

/*
 * interval binary search
 */

/*
inline int dsearch (double key, register double *v, int length) {

  register double result = 0;
  register int cursor;

  cursor = length /= 2;

  while (length != 0) {

    result = key - v [cursor];

    if      (result >   EPSILON) cursor += length/2;
    else if (result < - EPSILON) cursor -= length/2;
    else return cursor;

    if ((length>1) && (length % 2)) length++;
    length /= 2;
  }

  if (result < - EPSILON) return cursor - 1;
  else return cursor;
}
*/

 /*
  * block choice routine (vector linear)
  */

/*
void choose_block2 (sparseLP *lp,    // LP data
		   int *block,      // arrays with row indices to use (last is -1)
		   char *sat,       // fulfilled constraints
		   int nSat,        // number of fulfilled constraints
                   int cardinality, // cardinality of the block
		   double *b_Ax,    // violation vector
		   double sumviol   // sum of all violations
		   ) {

  register int i, j, rk = lp -> rk;

  double *random      = NULL;
  double  cum_rand    = 0;

  double *cum_weights = lp -> cum_weights;
  double  cum_weight  = 0;

  char   *violated    = lp -> chosen;

  //
  // generate cuWei
  //

  b_Ax += rk;

  for (j=0, i=rk; i>0; --i, --b_Ax)

    if (!(*sat++)) {
      cum_weight += exp (- *b_Ax / sumviol);
      *cum_weights++ = cum_weight;
      *violated++    = i;
      ++j;
    }

  cum_weights -= j;
  violated    -= j;
  sat         -= rk;

  //
  // generate y = a set of doubles (cumulatively)
  //

  for (i=cardinality; i>0; --i) {

    cum_rand += drand48 ();
    *random++ = cum_rand;
  }

  random -= cardinality;

  //
  // scale them down to cuWei [n]
  //

  for (i=cardinality; i>0; --i) {

    *random++ /= cum_rand;
  }

  random -= cardinality;

  // repeat 10 times
  //   for i=1..card
  //     linear search until right cuWei [k] found -> output k, card--


	int i, j;
	double cum_weight = 0.0;
	const int r = lp->rk;
	const int nViol = r - nSat;

	if (nViol == 0)
		return;

	for (i = 0, j = 0; i < r; ++i) {
		if (!sat[i]) {
			cum_weight += exp(-b_Ax[i]/sumviol);
		}
		cum_weights[i] = cum_weight;
	}

	for (i = 5 * cardinality - 1, j = 0; j < cardinality && i >= 0; --i) {

		const double key = drand48() * cum_weight;
		const int select = dsearch (key, cum_weights, r);
		if (chosen[select] == 0) {

			block[j] = select;
			chosen[select] = 1;
			++j;
		}
	}

	block [j--] = -1;

	for (; j >= 0; --j) {
		chosen[block[j]] = 0;
	}
}
*/
