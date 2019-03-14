// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
// Multiplies  two matrices. Matrices are stored in linear memory in row-major order,
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>


inline double dmin ( double a, double b ) { return a < b ? a : b; }
inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

// Matrix dimensions.


// Simple CPU implementation of matrix multiply
/* Parallelize with OpenMP */
  /* Pipeline and cache efficiency? */
// TO DO
/*     */

void MatrixMultiply(int rowsc, int colsc, int colsa, const double* A, const double* B, double* restrict C) 
{
  int rowc,colc,rowb,cola,idxa,idxb,idxc;
  for (rowc = 0; rowc < rowsc; ++rowc) {
    for (colc = 0; colc < colsc; ++colc) {
      double t = 0.0;
      for (cola=0; cola < colsa; ++cola) {
	idxa = rowc * colsa + cola;
	idxb = cola * colsc + colc;
	t += A[idxa] * B[idxb];
      }
      idxc = rowc * colsc + colc;
      C[idxc] = t;
    }
  }
}

const int ROWSC=100;
const int COLSA=100;
const int COLSC=100;


int main(int argc, char** argv) 
{

  /* Make sizes into arguments  from input */
  // TO DO
  /*     */


  
  // -----------------------  memory initialisation ----------------------- //
  
  double* A = (double*) malloc(sizeof(double)*ROWSC*COLSA);
  double* B = (double*) malloc(sizeof(double)*COLSA*COLSC);
  double* C = (double*) malloc(sizeof(double)*ROWSC*COLSC);
  int row, col, idx;

  /* Initialize with random data */
  // TO DO
  /*     */
  /* Transform to take NTRY measurements and keep the best */
  // TO DO 
  double t1 = omp_get_wtime();
  MatrixMultiply(ROWSC,COLSC,COLSA, A, B, C);
  double t2 = omp_get_wtime();
  double tmlt = (t2-t1);
  /*  */
  double mflops = (2.0e-6)*ROWSC*COLSC*COLSA/tmlt;
  fprintf(stdout,"Multiplying matrices of size %d x %d (%d) : time %lf  MFLOPS %lf \n",
	      ROWSC,COLSC,COLSA,tmlt,mflops);
  free(A);
  free(B);
  free(C);
  return 0;
}
