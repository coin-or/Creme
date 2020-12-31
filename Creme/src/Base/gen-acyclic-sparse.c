/* Builds a "sparse" acyclic subgraph instance: Parameters k,t,c
      Creates a k-inary tree with t+1 levels and replaces each node by a cycle with c nodes */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void write_sparse_problem (long k, long t, long c) {

  long i,j,l,parent,T,Tc,node,cyclevar;
  long Ei_minus, Ei,anz;

  if (c <= 1) {fprintf (stderr, "Error: c has to be larger than 1!"); exit(255);}

  T = (long) (pow (k, t+1) - 1) / (k-1);

  Tc = T * c;

  cyclevar = T - 1;

  //  fprintf ("Producing data set with:\n  %li variables\n  %li constraints\n  %li optimum\n",
  //	  Tc,Tc+T-1,Tc-1);

  // Write header and variables

  printf ("%li\t%li\t%li\t%li\n", Tc, Tc+T-1, (Tc+T-1) * 2, Tc-1);

  for (i=0; i<Tc; i++) printf ("0\t%li\n",Tc);

  anz      =  1;
  node     =  0;
  Ei_minus = -1;
  Ei       =  0;

  // insert cycle for node -1

  printf ("2\t1\t1e30\t%li\t1\t%i\t-1\n",cyclevar+1,0);

  for (l=0;l<c-2;l++) {
    cyclevar++;
    printf ("2\t1\t1e30\t%li\t1\t%li\t-1\n", cyclevar+1, cyclevar);
  }

  printf ("2\t1\t1e30\t%i\t1\t%li\t-1\n",0,cyclevar+1); cyclevar++;

  for (j=0; j<t; j++) {

    // One level of the tree

    for (i = anz*k; i>0; i--) {

      parent = Ei_minus + ((node - Ei) / k);

      printf ("2\t1\t1e30\t%li\t1\t%li\t-1\n", node+1, parent+1);  // write "forw-backw geq 1"

      // insert cycle:

      printf ("2\t1\t1e30\t%li\t1\t%li\t-1\n", cyclevar+1, node+1);

      for (l=0;l<c-2;l++) {
	cyclevar++;
	printf("2\t1\t1e30\t%li\t1\t%li\t-1\n",cyclevar+1,cyclevar);
      }

      printf ("2\t1\t1e30\t%li\t1\t%li\t-1\n",node+1,cyclevar+1);cyclevar++;

      ++node;
    }

    Ei_minus  = Ei;
    Ei        = node;
    anz      *= k;
  }
}

int main(int argc, char* argv[]){

  if (argc > 3) {

    const long k = atoi (argv [1]);
    const long t = atoi (argv [2]);
    const long c = atoi (argv [3]);

    write_sparse_problem (k, t, c);
  }
  else printf ("Builds a sparse acyclic subgraph instance\n\
Usage: %s k t c\n\
Creates a k-inary tree with t+1 levels and replaces each node by a cycle with c nodes\n", *argv);

  return 0;
}
