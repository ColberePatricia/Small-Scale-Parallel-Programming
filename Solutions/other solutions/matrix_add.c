// Derived from a CUDA example by Michał Czapiński 
// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Adds two matrices. Matrices are stored in linear memory in row-major order,
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include "wtime.h"

// Matrix dimensions.
const int ROWS = 4096;
const int COLS = 4096;


// Simple CPU implementation of matrix addition.
void MatrixAdd(int rows, int cols, const float* A, const float* B, float* C) 
{
  int row,col,idx;
#pragma omp parallel for private(row,col,idx) shared(A,B,C)
  for (row = 0; row < rows; ++row) {
    for (col = 0; col < cols; ++col) {
      idx = row * cols + col;
      C[idx] = A[idx] + B[idx];
    }
  }
}

int main(int argc, char** argv) 
{
  
  // ----------------------- Host memory initialisation ----------------------- //
  
  float* A = (float*) malloc(sizeof(float)*ROWS * COLS);
  float* B = (float*) malloc(sizeof(float)*ROWS * COLS);
  float* C = (float*) malloc(sizeof(float)*ROWS * COLS);
  int row, col, idx;
  
  srand(12345);
  for ( row = 0; row < ROWS; ++row) {
    for ( col = 0; col < COLS; ++col) {
      idx = row * COLS + col;
      A[idx] = 100.0f * ((float) rand()) / RAND_MAX;
      B[idx] = 100.0f * ((float) rand()) / RAND_MAX;
    }
  }
  
  double t1 = wtime();
  MatrixAdd(ROWS, COLS, A, B, C);
  double t2 = wtime();
  double tadd = (t2-t1);
  double mflops = (2.0e-6)*ROWS*COLS/tadd;
  fprintf(stdout,"Adding matrices of size %d x %d: time %lf  MFLOPS %lf \n",ROWS,COLS,tadd,mflops);
  free(A);
  free(B);
  free(C);
  return 0;
}
