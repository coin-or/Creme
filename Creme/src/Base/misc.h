/*
 * Name:   misc.h
 * Author: Pietro Belotti
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#ifndef MISC_H
#define MISC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <bzlib.h>

#include "sparse.h"
#include "lpio.h"
#include "rtr.h"

#define MAX_STR        40000
#define MALLOC_BLOCK    5000
#define CHUNKS_PER_LINE   10


void printLP (int *, sparseLP *); /* Print LP (for debugging purposes) */
void clearLP (sparseLP *);        /* Remove LP from memory             */

/*  Functions to read from a bzip compressed file */

int bzgetdbl  (BZFILE *, double *); /*  get double from .bz2 file */
int bzgetchar (BZFILE *, char *);   /*  get char   from .bz2 file */

static inline int bzgetint (BZFILE *b, register int *val) { /* get integer from .bz2 file */

  double x;
  register int ret_val;

  ret_val = bzgetdbl (b, &x);

  *val = floor (x);

  return ret_val;
}

/* Auxiliary functions to reallocate vectors of double/int */

static inline void reallocate_double (int req, int *size, double **buf) {

  if (req >= *size)
    *buf = (double *) realloc (*buf, (*size = (MALLOC_BLOCK * (req/MALLOC_BLOCK + 1))) * sizeof (double));
}

static inline void reallocate_int (int req, int *size, int **buf) {

  if (req >= *size)
    *buf = (int *) realloc (*buf, (*size = (MALLOC_BLOCK * (req/MALLOC_BLOCK + 1))) * sizeof (int));
}

void create_transpose (sparseLP *); /* Transpose coefficient matrix */


/* Get the norm of a vector of size l  */

static inline double get_norm (register double *v, register int l) {

  register double norm = 0;
  for (; l>0; --l, ++v) norm += *v * *v;
  return sqrt (norm);
}

#endif
