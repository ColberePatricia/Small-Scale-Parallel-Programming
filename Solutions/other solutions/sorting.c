#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"

inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

void MergeList(int *src1, int *src2, int len1, int len2, int *dest)
{
  int idx1=0, idx2=0, loc=0, i;
  while ( idx1 < len1 && idx2 < len2 )
    {
      if (src1[idx1] <= src2[idx2])
	{
	  dest[loc] = src1[idx1];
	  idx1++;
	}
      else
	{
	  dest[loc] = src2[idx2];
	  idx2++;
	}
      loc++;
    }

  for ( i=idx1; i<len1; i++)
    dest[loc++] = src1[i];

  for ( i=idx2; i<len2; i++)
    dest[loc++] = src2[i];
}

void MergeSort(int *data, int N)
{
  int *rep1, *rep2, *aux;

  int *temp=(int*)malloc(N*sizeof(int));
  int gSize, start, i;
  
  rep1=data;
  rep2=temp;

  for (gSize=1; gSize < N; gSize <<=1)
    {
      // Can you do better than this by using TASK? 
#pragma omp parallel for
      for ( start=0; start < N; start += 2*gSize)
	{
	  int next = start + gSize;
	  int g2Size = min(max(0,N-next),gSize);
	  if (g2Size == 0)
	    {
	      for ( i=0; i < N-start; i++)
		rep2[start+i] = rep1[start+i];
	    }
	  else
	    {
	      MergeList(rep1+start, rep1+next,gSize,g2Size,rep2+start);
	    }
	}
      aux=rep1; rep1=rep2; rep2=aux;
    }

  if (rep1 != data)
    memcpy(data,temp,N*sizeof(int));

  free(temp);
}


#define NUMS_PER_LINE  4
void prettyprint(int *data,int n, FILE *file)
{
  int i,j;
  for ( i=0; i<n; i+=NUMS_PER_LINE) {
    for ( j=0; j<NUMS_PER_LINE; j++) {
      if  (i+j < n) {
	fprintf(file,"%16d ",data[i+j]);
      }
    }
    fprintf(file,"\n");
  }
}

void main(int argc, char *argv[])
{
  if (argc <2) {
    fprintf(stderr,"sorting SIZE\n");
    exit(1);
  }
  int i;
  int N = atoi(argv[1]);
  int *data;

  if ((data=(int *) malloc(N*sizeof(int))) == NULL) {
    fprintf(stderr,"malloc failed\n");
    exit(1);
  }
  srand(12345);
  for ( i=0; i<N; i++)
    data[i] = rand();

  //fprintf(stdout,"Before sorting %d items: \n",N);
  //prettyprint(data,N,stdout);
  double t1=wtime();
  MergeSort(data,N);
  double t2=wtime();
  //fprintf(stdout,"After sorting %d items: \n",N);
  //prettyprint(data,N,stdout);
  fprintf(stdout,"Sorted %d items in %lf seconds\n",N,(t2-t1));
    
}
