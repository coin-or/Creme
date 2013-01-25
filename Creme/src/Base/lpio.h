/*
 * Name:   lpio.h
 * Author: Pietro Belotti
 *
 * This code is published under the Eclipse Public License (EPL).
 * See http://www.eclipse.org/legal/epl-v10.html
 *
 */

#ifndef LPIO_H
#define LPIO_H

#include "sparse.h"

int  read_problem (sparseLP *, char **);

void printLP      (int   *, sparseLP *);
void clearLP               (sparseLP *);

#endif
