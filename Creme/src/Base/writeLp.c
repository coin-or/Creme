/*
 * Name:    lpio.c
 * Author:  Pietro Belotti
 * Purpose: read compact .bz file into sparse matrix data structures
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <bzlib.h>

#include "sparse.h"
#include "lpio.h"
#include "rtr.h"
#include "misc.h"

/*
 * scan sequence of bzip2'ed files (each contains a subset of rows),
 * partition into chunks and send them
 */

int read_problem (sparseLP *lp,    /* LP data   */
 		  char **filename  /* file name */
		  ) {

  char first = 1;

  FILE *f;

  int status;

  BZFILE *bzf;

  int rownnz, curnnz, 
    totrow, currow = 0, 
    nnz, true_nnz, 
    i, j, best, 
    pos_i = 0, pos_d = 0, 
    size_d, size_i, 
    size_rl, size_ru,
    target_nnz,
    n_left_nnz;

  int bcast [4];
  int pos;

  static int    *bufi;
  static double *bufd;

  double lhsmax, lhsmin;

  double *rlb;
  double *rub;
  double *rhs;

  double  coe;

  double *lb, *ub;

  if (lp == NULL) {

    /*
      This function is called at the end with lp == NULL. The sole
      purpose of this call is to provide a mechanism to free bufi and
      bufd. That's why those variables are declared static, too.
    */

    if (bufi) free (bufi);
    if (bufd) free (bufd);

    return 0;
  }

  /*
   * Read chunks of matrix from .bz2 file 
   */

  if (lp -> my_id == 0) {

    lp->ntaut = lp->niis = 0;

    /*
     * init bz2'ed archive
     */

    if (strcmp (*filename, "-")) f = fopen (*filename, "r");
    else                         f = stdin;

    if ((f == NULL) || ferror (f)) {
      printf ("unable to open file %s\n", *filename);
      return -1;
    }
    
    printf ("%6.2f Opening bz2 file %s\n", CoinCpuTime(), *filename); fflush (stdout);

    bzf = BZ2_bzReadOpen (&status, f, 1, 0, NULL, 0);

    /*
     * read scalar parameters (file header)
     */
      
    bzgetint (bzf, &(lp -> c0));
    bzgetint (bzf, &(lp -> r0));
    bzgetint (bzf, &nnz);
    bzgetint (bzf, &best);

    n_left_nnz = true_nnz = nnz;

    /*
     * Send this initial information to slaves so that they prepare
     * the tables of the right size
     */

    bcast [0] = nnz;
    bcast [1] = lp -> c0;
    bcast [2] = lp -> r0;
    bcast [3] = lp -> ncpus;


    /*
     * create double and int buffers for the data to be sent to slaves
     */

    /*    target_nnz = nnz / lp -> ncpus + 1; */

    size_i = size_d = 0;
      
    bufd = NULL;
    bufi = NULL;

    /*
     * store variable upper/lower bounds (not alternate, but in
     * separate areas of the array)
     */

    for (i=lp->c0; i>0; i--, pos_d++) {

      reallocate_double (pos_d + lp -> c0, &size_d, &bufd);

      bzgetdbl (bzf, bufd + pos_d);
      bzgetdbl (bzf, bufd + pos_d + lp -> c0);
    }

    lb = bufd;
    ub = bufd + lp->c0;

    /*
     * allocate initial space for rhs (can be expanded with realloc)
     */

    size_rl = size_ru = 0;

    rlb = rub = NULL;

    /*
     * read row parameters and send chunk to each slave
     */

    /*    printf ("%6.2f: spreading matrix chunks...\n", CoinCpuTime()); */

    target_nnz = nnz / lp -> ncpus + 1;

    for (totrow = i = 0; i < lp -> ncpus; i++) {

      if (!(i % CHUNKS_PER_LINE)) {
	printf ("%6.2f chunks %d..%d\n", CoinCpuTime (), i, mymin (lp->ncpus-1, i + CHUNKS_PER_LINE - 1));
	/*	fflush (stdout); */
      }
      
      pos_d = 2 * lp->c0;
      pos_i = 0;

      /*
       * read each row and enqueue coefficients and indices
       */

      for (curnnz = currow = 0; 

	   ((i == lp->ncpus-1) || 
	    (curnnz < target_nnz)) && 
	     (totrow < lp -> r0); 

	   currow++, totrow++) {

	/*
	 * read this row's nnz
	 */

	bzgetint (bzf, &rownnz);

	/*	  printf ("%d: %d\n", currow, rownnz); fflush (stdout); */

	/*	printf ("line %d has %d nnzs (to be drawn from %d) %d < %d/%d\n",currow, rownnz, n_left_nnz); */

	n_left_nnz -= rownnz;
	
	curnnz += rownnz;

	reallocate_int (pos_i, &size_i, &bufi);

	bufi [pos_i++] = rownnz;

	reallocate_double (currow, &size_rl, &rlb);
	reallocate_double (currow, &size_ru, &rub);

	/*
	 * read b and c in b <= Ax <= c
	 */

	bzgetdbl (bzf, rlb + currow);
	bzgetdbl (bzf, rub + currow);

	if ((rlb [currow] > -1e20) && (rub [currow] < 1e20) && first) {

	  first = 0;
	  printf ("warning: range constraints (first found at %d: [%f,%f]) not yet implemented\n", 
		  currow, rlb [currow], rub [currow]);
	}

	/*
	 * read coefficients and indices
	 */

	lhsmin = lhsmax = 0;

#ifdef RTR_USE_PRAGMAS
#pragma execution_frequency (very_high)
#endif

	for (j=rownnz; j>0; --j, pos_i++, pos_d++) {

	  bzgetint (bzf, &pos);
	  if ((bzgetdbl (bzf, &coe) == 1) && (*++filename)) {

	    BZ2_bzReadClose (&status, bzf);
    
	    fclose (f);

	    if (strcmp (*filename, "-")) f = fopen (*filename, "r");
	    else                         f = stdin;

	    if ((f == NULL) || ferror (f)) {
	      printf ("unable to open file %s\n", *filename);
	      return -1;
	    }

	    printf ("%6.2f Opening bz2 file %s\n", CoinCpuTime(), *filename); fflush (stdout);

	    bzf = BZ2_bzReadOpen( &status, f, 1, 0, NULL, 0);

	    /*	      printf ("acabou!\n"); */
	  }

	  /* #pragma disjoint (*bufd, *bufi, coe, pos, pos_d, pos_i, lhsmax, lhsmin) */

	  {
	    reallocate_double (pos_d, &size_d, &bufd);
	    reallocate_int    (pos_i, &size_i, &bufi);

	    lb = bufd;
	    ub = bufd + lp->c0;

	    bufd [pos_d] = coe;
	    bufi [pos_i] = pos;

	    if (!(lp->noprep)) {

	      if (coe > 1e-20) {
		lhsmax += coe * ub [pos];
		lhsmin += coe * lb [pos];
	      }
	      else if (coe < -1e-20) {
		lhsmax += coe * lb [pos];
		lhsmin += coe * ub [pos];
	      }
	    }
	  }
	}

	if (!(lp->noprep)) {

	  if      (((rlb [currow] > -1e20) && (lhsmax <  rlb [currow])) || 
		   ((rub [currow] <  1e20) && (lhsmin >  rub [currow]))) { /* iis of size 1 */
	      
	    /*  printf ("iis %d: %.3e < %.3e or %.3e > %.3e\n", currow, lhsmax, rlb [currow], lhsmin, rub [currow]); */

	    pos_i -= (rownnz + 1); 
	    pos_d -=  rownnz; 
	    curnnz -= rownnz;
	    --currow;

	    ++ (lp -> niis);

	    true_nnz -= rownnz;
	    target_nnz = (true_nnz) / lp -> ncpus;
	  }
	  else if (((rub [currow] >  1e20) || (lhsmax <= rub [currow])) && 
		   ((rlb [currow] < -1e20) || (lhsmin >= rlb [currow]))) { /* tautology */

	    pos_i -= (rownnz + 1); 
	    pos_d -=  rownnz; 
	    curnnz -= rownnz;
	    --currow;

	    ++ (lp -> ntaut);
	    true_nnz -= rownnz;
	    target_nnz = (true_nnz) / lp -> ncpus;
	  }
	}
      }

      /*
       * Add constraints rhs's to the end of bufd
       */

      /* #pragma disjoint (*bufd, *rlb, *rub) */

      {
	reallocate_double (pos_d + 2 * currow, &size_d, &bufd);
	
	for (j=currow; j>0; j--) bufd [pos_d++] = *rlb++;
	for (j=currow; j>0; j--) bufd [pos_d++] = *rub++;
	
	rlb -= currow;
	rub -= currow;
      }

      lb = bufd;
      ub = bufd + lp->c0;

      bcast [0] = pos_d;
      bcast [1] = pos_i;
      bcast [2] = currow;

      bufd  = (double *) realloc (bufd, (size_d  = pos_d)  * sizeof (double));
      bufi  = (int    *) realloc (bufi, (size_i  = pos_i)  * sizeof (int));
      rlb   = (double *) realloc (rlb,  (size_rl = currow) * sizeof (double));
      rub   = (double *) realloc (rub,  (size_ru = currow) * sizeof (double));

      lb = bufd;
      ub = bufd + lp->c0;

      /*
      if      (currow < 1e4)  printf ("[%d,",  currow);
      else if (currow < 1e7)  printf ("[%dk,", currow/(int)1e3);
      else if (currow < 1e10) printf ("[%dm,", currow/(int)1e6);
      else                    printf ("[%dg,", currow/(int)1e9);

      if      (pos_i <= 1e4) printf ("%d] ",   pos_i - currow);
      else if (pos_i <= 1e7) printf ("%dk] ", (pos_i - currow)/(int)1e3);
      else if (pos_i <= 1e7) printf ("%dm] ", (pos_i - currow)/(int)1e6);
      else                   printf ("%dg] ", (pos_i - currow)/(int)1e9); fflush (stdout);

      if (!((i+1) % CHUNKS_PER_LINE)) printf ("\n");
      */

      if (i < lp->ncpus - 1) 
	target_nnz = n_left_nnz / (lp->ncpus - i - 1);
    }

    printf   ("%6.2f done. %d rows, %d columns, %d nonzero", CoinCpuTime (), lp->r0, lp->c0, nnz);
      
    if (!(lp -> noprep) && (lp->niis || lp->ntaut)) 
      printf ("\n               %d iis, %d tautologies, %d nonzero eliminated", lp->niis, lp->ntaut, nnz - true_nnz);

    printf ("\n");

    free (rlb);
    free (rub);

    BZ2_bzReadClose (&status, bzf);

    fclose (f);
  }

  /*
   * create sparse representation of local submatrix A_j
   */

  lp -> rk = currow;
  lp -> nnzk = pos_d - 2 * (lp -> rk + lp -> c0);

  lp -> lb  = bufd;
  lp -> ub  = lp -> lb + lp -> c0;

  lp -> rlb = lp -> ub + lp -> c0 + lp -> nnzk;
  lp -> rub = lp -> rlb + lp -> rk;

  rhs = lp -> rhs = (double *) malloc (currow * sizeof (double));

  rlb = lp -> rlb;
  rub = lp -> rub;

  for (j=currow; j>0; j--, rub++, rlb++, rhs++) 

    if (*rlb < -1e20) *rhs = *rub;
    else              *rhs = *rlb;

  rhs -= currow;
  rlb -= currow;
  rub -= currow;

  lp -> ic = (double **) malloc (currow * sizeof (double *));
  lp -> ip = (int    **) malloc (currow * sizeof (int    *));
  lp -> il = (int     *) malloc (currow * sizeof (int));

  lp -> chosen = (char *) malloc (currow * sizeof (char));

  pos_i = 0;
  pos_d = 2 * lp -> c0;

  {
    double **ic = lp -> ic;
    int    **ip = lp -> ip;
    int     *il = lp -> il;
    char    *ch = lp -> chosen;

    for (j=currow; j>0; j--) {

      *ch++ = 0;
      *ic++ = bufd + pos_d;
      *ip++ = bufi + pos_i + 1;
      *il++ = bufi [pos_i];

      pos_d += bufi [pos_i];
      pos_i += bufi [pos_i] + 1;
    }

    ic -= currow;
    il -= currow;

    /*
     * Normalization
     */

    for (i=lp->rk; i>0; i--, ic++, il++, rlb++, rub++, rhs++) {

      register double norm = get_norm (*ic, *il);

      if (norm > EPSILON) {

	register double *icl = *ic;
      
	if (*rub < 1e29) {

	  norm = - norm;
	  *rhs = *rub / norm;
	}
	else *rhs = *rlb / norm;

	for (j = *il; j > 0; j--) *icl++ /= norm;
      }
    }
  }

  /*  if (lp->my_id==0) printf ("%6.2f create transposed matrix\n", CoinCpuTime ()); */

  create_transpose (lp);

  return status;
}
