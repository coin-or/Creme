/*
 * Name:    rtr.h
 * Author:  Pietro Belotti
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 */

#ifndef RTR_H
#define RTR_H

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include "sparse.h"

enum {USE_RTR = 0, USE_LOCSRCH};

/*
 * input: the unfeasible LP
 * output: a char vector containing a MFS
 */

int rtr (sparseLP *lp, char *sol);

/*
 * Returns cpu time
 */

static inline double CoinCpuTime()
{
  double cpu_temp;

#if defined(_MSC_VER)
  unsigned int ticksnow;        /* clock_t is same as int */
  
  ticksnow = (unsigned int)clock();
  
  cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  cpu_temp = usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
#endif

  return cpu_temp;
}

#endif
