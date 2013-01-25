/*
 * Name:    isfeas.h
 * Author:  Pietro Belotti
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#ifndef ISFEAS_H
#define ISFEAS_H

#include "sparse.h"

int isFeas (sparseLP *,    /* sparse LP data                          */
	    char     *,    /* alleged fulfilled constraints           */
	    double   *,    /* current point                           */
            int      *);   /* alleged number of constraints satisfied */

#endif
