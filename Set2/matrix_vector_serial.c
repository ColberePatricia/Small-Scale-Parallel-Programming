// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>


inline double dmin ( double a, double b ) { return a < b ? a : b; }
inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

// Matrix dimensions.
const int ROWS = 4096;
const int COLS = 4096;


// Simple CPU implementation of matrix-vector product
void MatrixVector(int rows, int cols, const double* A, const double* x, double* restrict y) 
{
  /* Parallelize with OpenMP */
  /* Pipeline and cache efficiency? */
  // TO DO
  /*     */
  int row,col, idx;
  double t;
  for (row = 0; row < rows; ++row) {
    t=0.0;
    for (col = 0; col < cols; ++col) {
      idx = row * cols + col;
      t = t + A[idx]*x[col];
    }
    y[row] = t;
  }
}

int main(int argc, char** argv) 
{
  /* Make ROWS and COLS into arguments  from input */
  // TO DO
  /*     */

  
  double* A = (double*) malloc(sizeof(double)*ROWS * COLS);
  double* x = (double*) malloc(sizeof(double)*ROWS );
  double* y = (double*) malloc(sizeof(double)*ROWS );
  int row, col, idx;
  /* Initialize with random data */
  // TO DO
  /*     */
    

  /* Transform to take NTRY measurements and keep the best */
  // TO DO 
  double t1 = omp_get_wtime();
  MatrixVector(ROWS, COLS, A, x, y);
  double t2 = omp_get_wtime();
  double tmlt = (t2-t1);
  /*   */ 
  double mflops = (2.0e-6)*ROWS*COLS/tmlt;
  
  fprintf(stdout,"Matrix-Vector product of size %d x %d with 1 thread: time %lf  MFLOPS %lf \n",
	  ROWS,COLS,tmlt,mflops);
  free(A);
  free(x);
  free(y);
  return 0;
}
