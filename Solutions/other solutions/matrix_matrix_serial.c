// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
// Multiplies  two matrices. Matrices are stored in linear memory in row-major order,
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"
const int ntimes = 5;
inline double dmin ( double a, double b ) { return a < b ? a : b; }
inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

// Matrix dimensions.


// Simple CPU implementation of matrix multiplication.
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

int main(int argc, char** argv) 
{
  if (argc < 4) {
    fprintf(stderr,"Usage: %s  rowsc colsc colsa\n",argv[0]);
  }
  int nrowsc=atoi(argv[1]);
  int ncolsc=atoi(argv[2]);
  int ncolsa=atoi(argv[3]);
  
  
  // ----------------------- Host memory initialisation ----------------------- //
  
  double* A = (double*) malloc(sizeof(double)*nrowsc*ncolsa);
  double* B = (double*) malloc(sizeof(double)*ncolsa*ncolsc);
  double* C = (double*) malloc(sizeof(double)*nrowsc*ncolsc);
  int row, col, idx;
  
  srand(12345);
  for ( row = 0; row < nrowsc; ++row) {
    for ( col = 0; col < ncolsa; ++col) {
      idx = row * ncolsa + col;
      A[idx] = 100.0f * ((double) rand()) / RAND_MAX;
    }
  }
  for ( row = 0; row < ncolsa; ++row) {
    for ( col = 0; col < ncolsc; ++col) {
      idx = row * ncolsc + col;
      B[idx] = 100.0f * ((double) rand()) / RAND_MAX;
    }
  }

  double tmlt = 1e100;
  for (int try=0; try < ntimes; try ++ ) {
    double t1 = wtime();
    MatrixMultiply(nrowsc,ncolsc,ncolsa, A, B, C);
    double t2 = wtime();
    tmlt = dmin(tmlt,(t2-t1));
  }
  double mflops = (2.0e-6)*nrowsc*ncolsc*ncolsa/tmlt;
  fprintf(stdout,"Multiplying matrices of size %d x %d (%d) : time %lf  MFLOPS %lf \n",
	      nrowsc,ncolsc,ncolsa,tmlt,mflops);
  free(A);
  free(B);
  free(C);
  return 0;
}
