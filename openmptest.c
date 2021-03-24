#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

/* A simple example of using a few of the main constructs in openmp.
 * To compile this (on linux):
 *   gcc -o openmptest -fopenmp openmptest.c
 */
int main(int argc, char *argv[]) {
  int i;
  int t = 8, n = 10;
  if (argc > 1) {
    n = atoi(argv[1]);
    printf("Using %d iterations\n", n);
    if (argc > 2) {
      t = atoi(argv[2]);
      printf("Using %d threads\n", t);
    }
  }
  /* Use the openmp api to set the max number of threads dynamically for
     the next parallel region. */
  omp_set_dynamic(0); /* Disable dynamic teams. */
  omp_set_num_threads(t);

  /* A parallel for loop, iterations divided amongst the threads */
#pragma omp parallel for
  for (i = 0; i < n; i++) {
    printf("Iteration %d done by thread %d\n",
           i,
           omp_get_thread_num());
  }

  /* A parallel stmt, executed by all threads */
  /* Note that we didn't specify the number of threads here. */
#pragma omp parallel
  printf("Hello from thread %d, nthreads %d\n",
         omp_get_thread_num(),
         omp_get_num_threads());

  /* Three parallel sections doing different work */
#pragma omp parallel
#pragma omp sections
  {
#pragma omp section
    {
      int j, k = 0;
      for (j = 0; j < 100000; j++)
        k += j;
      printf("Section 1 from thread %d, nthreads %d, val=%d\n",
             omp_get_thread_num(),
             omp_get_num_threads(), k);
    }
#pragma omp section
    {
      int j, k = 0;
      for (j = 0; j < 100000; j++)
        k += j;
      printf("Section 2 from thread %d, nthreads %d, val=%d\n",
             omp_get_thread_num(),
             omp_get_num_threads(), k);
    }
#pragma omp section
    {
      int j, k = 0;
      for (j = 0; j < 100000; j++)
        k += j;
      printf("Section 3 from thread %d, nthreads %d, val=%d\n",
             omp_get_thread_num(),
             omp_get_num_threads(), k);
    }
  }
}
