// Derived from a CUDA example by Michał Czapiński 
// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Adds two matrices. Matrices are stored in linear memory in row-major order,
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

// Matrix dimensions.
const int ROWS = 4096;
const int COLS = 4096;


// Simple CPU implementation of matrix addition.
void MatrixAdd(int rows, int cols, const float* A, const float* B, float* restrict C) 
{
  /* Parallelize with OpenMP */
  // TO DO
  /*     */
  int row,col,idx;
  for (row = 0; row < rows; ++row) {
    for (col = 0; col < cols; ++col) {
      idx = row * cols + col;
      C[idx] = A[idx] + B[idx];
    }
  }
}

int main(int argc, char** argv) 
{
  /* Make ROWS and COLS into arguments  from input */
  // TO DO
  /*     */

  
  // -----------------------  memory initialisation ----------------------- //
  
  float* A = (float*) malloc(sizeof(float)*ROWS * COLS);
  float* B = (float*) malloc(sizeof(float)*ROWS * COLS);
  float* C = (float*) malloc(sizeof(float)*ROWS * COLS);
  int row, col, idx;
  

  /* Initialize with random data */
  // TO DO
  /*     */

  double t1 = omp_get_wtime();
  MatrixAdd(ROWS, COLS, A, B, C);
  double t2 = omp_get_wtime();
  double tadd = (t2-t1);
  double mflops = (1.0e-6)*ROWS*COLS/tadd;
  fprintf(stdout,"Adding matrices of size %d x %d: time %lf  MFLOPS %lf \n",ROWS,COLS,tadd,mflops);
  free(A);
  free(B);
  free(C);
  return 0;
}
