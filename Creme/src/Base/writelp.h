/*
 * Name:    writelp.c
 * Author:  Pietro Belotti
 * Purpose: write feasible subsystem into a .lp file
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#ifndef WRITELP_H
#define WRITELP_H

#include "sparse.h"

/*
 * Write feasible subsystem into a .lp file
 */

int write_lp (sparseLP *lp,   /* LP data   */
	      char *mfs,      /* feasible subsystem */
	      char *filename  /* file name */
	      );

#endif
