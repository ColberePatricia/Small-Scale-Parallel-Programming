#include <omp.h>
#include <stdio.h>
#define N  10000
main()
{
  int i,k,n,tid;
  float x[N], y[N];

  float dot=0.0;
  n=N;
  
  for (i=0; i<n; i++) {
    x[i] = 2.0*i;
    y[i] = 3.0*i;
  }
#if 0
  /* this is functionally correct */
#pragma omp parallel shared(n,x,y,dot) private(tid)
  {
    tid = omp_get_thread_num();
    float mydot = 0.0;
#pragma omp for
    for (i=0; i<n; i++)
      mydot += x[i]*y[i];
#pragma omp critical (update_dot)
    {
      dot += mydot;
      printf("thread %d: mydot %f  dot %f\n",
             tid,mydot,dot);      
    }
  }
#else
    /* this is the best solution, the compiler
       will optimize the summation pattern */
#pragma omp parallel for shared(n,x,y) reduction(+:dot)
    for (i=0; i<n; i++)
      dot += x[i]*y[i];
#endif
  printf("size %d   dot %f\n",n,dot);
}
