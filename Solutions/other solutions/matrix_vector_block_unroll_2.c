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

void innerMatrixVector(int rows, int cols, const double* A, int ncallc, const double* x, 
		       double beta, double* restrict y)
{
  int row,col, idx;
  double t0, t1;

#pragma omp parallel for shared(x,y,A,rows,cols,ncallc) private(row,col,idx,t0, t1 ) 
  for (row = 0; row < rows - rows%2; row += 2) {
    if (beta == 0.0) {
      t0 = beta;
      t1 = beta;
    } else {
      t0 = beta *  y[row+0];
      t1 = beta *  y[row+1];
    }

    for (col = 0; col < cols - cols%4 ; col+=4) {
      t0 += A[(row+0) * ncallc + col+0] * x[col+0] + A[(row+0) * ncallc + col+1] * x[col+1]
	   +A[(row+0) * ncallc + col+2] * x[col+2] + A[(row+0) * ncallc + col+3] * x[col+3];
      t1 += A[(row+1) * ncallc + col+0] * x[col+0] + A[(row+1) * ncallc + col+1] * x[col+1]
	   +A[(row+1) * ncallc + col+2] * x[col+2] + A[(row+1) * ncallc + col+3] * x[col+3];
    }
    for (col = cols - cols%4; col < cols; col++) {
      t0 += A[(row+0) * ncallc + col] * x[col] ;
      t1 += A[(row+1) * ncallc + col] * x[col] ;
    }
    
    y[row+0] = t0;
    y[row+1] = t1;
  }
  
  for (row = rows - rows%2; row < rows; row++) {
    double t=0.0;
    for (col = 0; col < cols; col++) {
      int idx = row * ncallc + col;
      t = t + A[idx]*x[col];      
    }
    y[row]=t;
  }
}



const int BBS=1000;
void MatrixVector(int rows, int cols, const double* A, const double* x, double* restrict y) 
{
  int row,col;
  
  for (row = 0; row < rows; row += BBS) {
    int nr  = min(BBS, rows-row);
    col     = 0;
    int nc  = cols;
    int idx = row * cols + col;    
    innerMatrixVector(nr, nc, &(A[idx]), cols, &(x[col]), 0.0, &(y[row]));
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
      fprintf(stdout,"Matrix-Vector product (block_unroll_2) of size %d x %d with %d threads: time %lf  MFLOPS %lf \n",
	      nrows,ncols,omp_get_num_threads(),tmlt,mflops);
    }
  }
  free(A);
  free(x);
  free(y);
  return 0;
}
