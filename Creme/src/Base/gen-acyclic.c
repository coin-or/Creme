#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void write_acyclic_problem(long n_nodes, long n_add_edges) {

  long i,j, i2,j2, k, n,N,
       n_forw = n_nodes * (n_nodes-1) / 2;

  if (n_add_edges) N = n_forw / n_add_edges;
  else             N = n_forw + 1;

  printf ("%li\t%li\t%li\t%li\n",
	  n_nodes,
	  n_forw + n_add_edges,
	  (n_forw+n_add_edges) * 2,
	  n_forw);

  //
  // variables
  //

  for (i=0;i<n_nodes;i++) printf("0\t%li\n",n_nodes);

  //
  // constraints
  //

  for (n=k=i=0; i<n_nodes; i++) {

    for (j=i+1; j<n_nodes; j++) {

      if (!(n++ % N) && (k++ < n_add_edges)) {

	//
	// add additional backward arcs
	//

	do {

	  i2 = floor (drand48 () * n_nodes);
	  j2 = floor (drand48 () * n_nodes);

	} while (i2==j2);

	if (i2 < j2) printf("2\t1\t1e30\t%li\t1\t%li\t-1\n", j2, i2);
	else         printf("2\t1\t1e30\t%li\t1\t%li\t-1\n", i2, j2);
      }

      printf("2\t1\t1e30\t%li\t1\t%li\t-1\n",i,j);
      // forw-backw geq 1, i.e.,i-j
    }
  }
}

int main(int argc, char* argv[]){

  if (argc >= 2) {

    const int n_nodes     = atoi (argv [1]);
    const int n_add_edges = (argc >= 3) ? (atoi (argv [2])):0;

    write_acyclic_problem (n_nodes, n_add_edges);
  }
  else printf ("usage: %s <num_nodes> [<num_backward_arcs>]\n", *argv);

  return 0;
}
