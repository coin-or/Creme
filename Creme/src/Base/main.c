/*
 * Name:    main.c
 * Author:  Pietro Belotti
 * Purpose: define program options, read problem, call rtr procedure
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <signal.h>
#include <string.h>

#include "sparse.h"
#include "lpio.h"
#include "writelp.h"
#include "rtr.h"
#include "cmdline.h"

int interrupt = 0;

void userInterrupt () {

  if (interrupt++ > 1) {
    fprintf (stderr, "User interrupt, aborting.\n");
    exit (-1);
  }
}

/*
 * Main
 */

int main (int argc, char **argv) {

  char **filenames;

  sparseLP lp;

  char *mfs, help = 0;

  char *output = (char *) malloc (1 * sizeof (char));

  struct timeval tv;

  int i;
  int rndseed = 1; /* random seed fixed to 1 if user doesn't set it */

  /*
   * Specify program command line options
   */

  tpar options []   = {{'1', "1dim",       0,      NULL,    TTOGGLE, "use 1-dimensional optimization"}
		      ,{'f', "noprep",     0,      NULL,    TTOGGLE, "no pre-processing on vars/cons"}
		      ,{'M', "bigm",       0,      NULL,    TTOGGLE, "just create the big-M MILP and exit"}
		      ,{'U', "normalize",  0,      NULL,    TTOGGLE, "normalize best point found"}
		      ,{'i', "inv-cool",   0,      NULL,    TTOGGLE, "inverse linear temperature decrease"}
		      ,{'l', "lin-cool",   0,      NULL,    TTOGGLE, "linear temperature decrease"}
		      ,{'d', "dbl-rand",   0,      NULL,    TTOGGLE, "use double randomization"}

		      ,{'a', "alpha",      2,      NULL,    TDOUBLE, "set temperature scaling factor"}
		      ,{'b', "beta",       0.5,    NULL,    TDOUBLE, "set convex combination factor"}
		      ,{'g', "gammarate",  1.0003, NULL,    TDOUBLE, "set temperature decreasing factor"}
		      ,{'m', "murate",     1.02,   NULL,    TDOUBLE, "set block shrinking in worsening iters."}
		      ,{'e', "stretch",    1e4,    NULL,    TDOUBLE, "set stretch on d [t] in 1-dim. opt."}

		      ,{'n', "iter",       100000, NULL,    TINT,    "set max. # iterations"}
		      ,{'s', "seed",       0,      NULL,    TINT,    "set random seed"}

		      ,{'t', "time",       -1,     NULL,    TDOUBLE, "set maximum runtime"}
		      ,{'r', "restart",    0,      NULL,    TINT,    "restart after <num> iterations"}
		      ,{'D', "dualsearch", 0,      NULL,    TINT,    "variable local search on <n> coords"}
		      ,{'i', "infinity",   1e40,   NULL,    TDOUBLE, "value for infinity"}

	              ,{'o', "output",     0,      NULL,    TSTRING, "filename for feasible subsystem"}
		      ,{'h', "help",       0,      NULL,    TTOGGLE, "print this help"}

  	              ,{0, "", 0, NULL, 0, ""}            /* THIS ENTRY ALWAYS AT THE END */
  };

  options  [0].par =  &(lp.onedim);
  options  [1].par =  &(lp.noprep);
  options  [2].par =  &(lp.bigm);
  options  [3].par =  &(lp.norm);
  options  [4].par =  &(lp.invcool);
  options  [5].par =  &(lp.lincool);
  options  [6].par =  &(lp.dblrand);

  options  [7].par =  &(lp.alpha);
  options  [8].par =  &(lp.beta);
  options  [9].par =  &(lp.gammaRate);
  options [10].par =  &(lp.muRate);
  options [11].par =  &(lp.stretch);

  options [12].par =  &(lp.nIter);
  options [13].par =  &rndseed;

  options [14].par =  &(lp.timelimit);
  options [15].par =  &(lp.restFreq);
  options [16].par =  &(lp.blkcard);
  options [17].par =  &(lp.infinity);

  options [18].par =  &output;
  options [19].par =  &help;

  options [20].par =  NULL;

  /*
   * default parameter values
   */

  set_default_args (options);

  lp.my_id = 0;
  lp.ncpus = 1;

  /*
   * parse command line
   */

  /*  if (lp.my_id == 0) printf ("reading input\n"); */

  filenames = readargs (argc, argv, options);

  /*
   * initialize random seed (modifiable through command line)
   */

  if (rndseed == 0) { /* do it like AMPL */

    gettimeofday (&tv, NULL);
    rndseed = tv.tv_sec * 1e6 + tv.tv_usec;
    printf ("Random seed initialized at %d\n", rndseed);
  }

  if (help) print_help (*argv, options);

  if (!filenames) {

    if (!help) {

      printf ("%s: no input given\n", *argv);
      print_help (*argv, options);
      return -1;
    }
    else return 0;
  }

  lp.locsea = (lp.blkcard != 0);

  /*
   * random generator startup
   */

  srand48 (rndseed + lp.my_id);

  signal (SIGINT, userInterrupt);

  /*
   * read instance
   */

  if (read_problem (&lp, filenames) < 0)
    return -1;

  mfs = (char *) malloc (lp.rk * sizeof (char));

  if (lp.restFreq <= 0) lp.restFreq = lp.nIter;

  /*
   * run the R.T.R.
   */

  if (lp.my_id == 0) printf ("%6.2f start rtr\n", CoinCpuTime());

  if (rtr (&lp, mfs)) {

    /*    if (lp.lpSolv) solveFeasLP   (olp, &lp, mfs); */
    if (strcmp (output, "")) write_lp (&lp, mfs, output);
  }
  else printf ("Error in solving the problem\n");

  read_problem (NULL, NULL);

  for (i=0; filenames [i] != NULL; i++) free (filenames [i]);

  free (filenames);
  free (mfs);
  free (output);

  clearLP (&lp);

  return 0;
}
