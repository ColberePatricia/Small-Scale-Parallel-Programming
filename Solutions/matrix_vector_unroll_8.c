// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"
const int ntimes=50;
inline double dmin ( double a, double b ) { return a < b ? a : b; }

inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

// Unrolling by 8 rows
void MatrixVector(int rows, int cols, const double* A, const double* x, double* restrict y) 
{
  int row,col, idx;
#pragma omp parallel for shared(x,y,A) private(row,col,idx )
  for (row = 0; row < rows - rows%8; row += 8) {
    double t0 = A[(row+0) * cols] * x[0];
    double t1 = A[(row+1) * cols] * x[0];
    double t2 = A[(row+2) * cols] * x[0];
    double t3 = A[(row+3) * cols] * x[0];
    double t4 = A[(row+4) * cols] * x[0];
    double t5 = A[(row+5) * cols] * x[0];
    double t6 = A[(row+6) * cols] * x[0];
    double t7 = A[(row+7) * cols] * x[0];
    for (col = 1; col < cols ; col++) {
      t0 += A[(row+0) * cols + col] * x[col];
      t1 += A[(row+1) * cols + col] * x[col];
      t2 += A[(row+2) * cols + col] * x[col];
      t3 += A[(row+3) * cols + col] * x[col];
      t4 += A[(row+4) * cols + col] * x[col];
      t5 += A[(row+5) * cols + col] * x[col];
      t6 += A[(row+6) * cols + col] * x[col];
      t7 += A[(row+7) * cols + col] * x[col];
    }
    y[row+0] = t0;
    y[row+1] = t1;
    y[row+2] = t2;
    y[row+3] = t3;
    y[row+4] = t4;
    y[row+5] = t5;
    y[row+6] = t6;
    y[row+7] = t7;
  }
  
  for (row = rows - rows%8; row < rows; row++) {
    double t=0.0;
    for (col = 0; col < cols; col++) {
      int idx = row * cols + col;
      t = t + A[idx]*x[col];      
    }
    y[row]=t;
  }
}


int main(int argc, char** argv) 
{
  
  if (argc < 3) {
    fprintf(stderr,"Usage: %s  rows cols\n",argv[0]);
  }
  int nrows=atoi(argv[1]);
  int ncols=atoi(argv[2]);
  


  
  double* A = (double*) malloc(sizeof(double)*nrows * ncols);
  double* x = (double*) malloc(sizeof(double)*nrows );
  double* y = (double*) malloc(sizeof(double)*nrows );
  int row, col, idx;
  
  srand(12345);
  for ( row = 0; row < nrows; ++row) {
    for ( col = 0; col < ncols; ++col) {
      idx = row * ncols + col;
      A[idx] = 100.0f * ((double) rand()) / RAND_MAX;
    }
    x[row] = 100.0f * ((double) rand()) / RAND_MAX;      
  }
  
  double tmlt = 1e100;
  for (int try=0; try < ntimes; try ++ ) {
    double t1 = wtime();
    MatrixVector(nrows, ncols, A, x, y);
    double t2 = wtime();
    tmlt = dmin(tmlt,(t2-t1));
  }
  double mflops = (2.0e-6)*nrows*ncols/tmlt;
#pragma omp parallel 
  {
#pragma omp master
    {
      fprintf(stdout,"Matrix-Vector product (unroll_8) of size %d x %d with %d threads: time %lf  MFLOPS %lf \n",
	      nrows,ncols,omp_get_num_threads(),tmlt,mflops);
    }
  }
  free(A);
  free(x);
  free(y);
  return 0;
}
