#include "wtime.h"

double wtime() 
{
  clock_t t=clock();
  return(   ((double) t/CLOCKS_PER_SEC) );
}
