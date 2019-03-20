// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include "wtime.h"
#include "mmio.h"
#include "matrixPreprocessing.h"
#include "test.h"


// Generates a random vector
void generateVector(int vectorSize, double* vector) {
	srand(12345);
	for (int row = 0; row < vectorSize; ++row) {
		vector[row] = 100.0f * ((double)rand()) / RAND_MAX;
	}
}

// Generates a random matrix
void generateMatrix(int rows, int cols, double* matrix) {
	srand(21345);
	int idx;
	for (int row = 0; row < rows; ++row) {
		for (int col = 0; col < cols; ++col) {
			idx = row * cols + col;
			matrix[idx] = 100.0f * ((double)rand()) / RAND_MAX;
		}
	}
}


// Simple CPU implementation of matrix-vector product
void MatrixVector(int rows, int cols, const double* A, const double* x, double* y) 
{
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

void MatrixVectorCSR(int M, int N, double* IRP, double* JA, double* AS, double* x, double* y)
{
	/*
	for i=1:m
	t=0;
	for j=irp(i):irp(i+1)-1
	t = t + as(j)*x(ja(j));
	2
	end
	y(i) = t;
	end
	*/
}

void MatrixVectorELLPACK(int M, int N, int MAXNZ, double* JA, double* AS, double* x, double* y)
{
	/*
	for i=1:m
	t=0;
	for j=1:maxnzr
	t = t + as(i,j)*x(ja(i,j));
	end
	y(i) = t;
	end
	*/
}

int main(int argc, char** argv) 
{
  
	testMatrixProcessing();

	char* matrixFile;
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\nWe will be using the cage4 matrix as default\n", argv[0]);
		matrixFile = "D:\\Cranfield work\\Small Scale Parallel Programming\\matrices\\cage4.mtx";
	}
	else {
		matrixFile = argv[1];
	}

	readMatrix(matrixFile);

	//------------------
	int nrows=2;
	int ncols=2;
	/*fprintf(stdout, "Number of rows: ");
	fscanf_s(stdin, "%d", &nrows);
	fprintf(stdout, "\nNumber of columns: ");
	fscanf_s(stdin, "%d", &ncols);  
	fprintf(stdout, "\n");*/


  
  double* A = (double*) malloc(sizeof(double)*nrows * ncols);
  double* x = (double*) malloc(sizeof(double)*nrows );
  double* y = (double*) malloc(sizeof(double)*nrows );
  
  generateVector(nrows, x);
  generateMatrix(nrows, ncols, A);

  /*
  // We check the matrix and the vector
  fprintf(stdout, "Matrix:\n");
  printMatrix(nrows, ncols, A);
  fprintf(stdout, "Vector:\n");
  printVector(nrows, x);
  */

  double t1 = wtime();
  MatrixVector(nrows, ncols, A, x, y);
  double t2 = wtime();
  double tmlt = (t2-t1);
  double mflops = (2.0e-6)*nrows*ncols/tmlt;
  
  fprintf(stdout,"Matrix-Vector product of size %d x %d with 1 thread: time %lf  MFLOPS %lf \n", nrows,ncols,tmlt,mflops);
  free(A);
  free(x);
  free(y);
  return 0;
}
