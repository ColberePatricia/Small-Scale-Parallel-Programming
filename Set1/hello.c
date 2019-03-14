#include <omp.h>
#include <stdio.h>

int main () {
  int nthreads, tid;
  /* Fork a team of threads with each 
     thread having a private tid variable */
  // TO DO
  /*   */ 

  /* Obtain and print thread id */
  // TO DO    
  /*   */ 

  /* Only master thread does this */
  if (tid == 0)  {
    /* obtain number of threads */
    // TO DO
    /*   */ 
    printf("Number of threads = %d\n", nthreads);
  }
  
  /* All threads join master thread and terminate */
  return 0;
}
