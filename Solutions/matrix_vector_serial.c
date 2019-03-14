// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include "wtime.h"
#include "mmio.h"


// Print the matrix in the console
void printMatrix(int rows, int cols, const double* A) {
	int idx;
	for (int row = 0; row < rows; ++row) {
		for (int col = 0; col < cols; ++col) {
			idx = row * cols + col;
			fprintf(stdout, "%lf, ", A[idx]);
		}
		fprintf(stdout, "\n");
	}
}

// Print the vector in the console
void printVector(int size, const double* vector) {
	for (int idx = 0; idx < size; ++idx) {
		fprintf(stdout, "%lf, ", vector[idx]);
	}
	fprintf(stdout, "\n");
}

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

void convertToNonCompressedMatrix(int M, int N, int nz, int *I, int *J, double* val, double* result) {
	// We reserve memory for the result matrix
	result = (double *)malloc(M * N * sizeof(double));
	// We initialize the matrix with zeroes
	for (int i = 0;i < M*N;i++)
		result[i] = 0;

	int idx;

	for (int i = 0;i < nz;i++) {
		idx = I[i]*N + J[i];
		//fprintf(stdout, "I=%d, J=%d, idx=%d\n", I[i], J[i], idx);
		result[idx] = val[i];
	}
}

void convertToELLPACK(int M, int N, int nz, int *I, int *J, double* val, double* result) {

}

void convertToCSR(int M, int N, int nz, int *I, int *J, double* val, double* result) {

}

// Read a matrix from the input
void readMatrix(char* fileName) {
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;
	int i, *I, *J;
	double *val;

	// If the file of the matrix cannot be opened
	if ((f = fopen(fileName, "r")) == NULL)
		fprintf(stdout, "The file %s could not be opened\n", fileName);
	else if (mm_read_banner(f, &matcode) != 0)
		fprintf(stdout, "Could not process Matrix Market banner.\n");
	else if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode))
		fprintf(stdout, "Sorry, this application does not support Market Market type: [%s]\n", mm_typecode_to_str(matcode));
	else if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
		fprintf(stdout, "Could not read the size of the matrix");
	else {
		/* reseve memory for matrices */
		I = (int *)malloc(nz * sizeof(int));
		J = (int *)malloc(nz * sizeof(int));
		val = (double *)malloc(nz * sizeof(double));


		/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
		/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
		/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
		for (i = 0; i < nz; i++)
		{
			fscanf_s(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
			I[i]--;  // adjust from 1-based to 0-based
			J[i]--;
		}

		if (f != stdin) fclose(f);

		/************************/
		/* now write out matrix */
		/************************/

		mm_write_banner(stdout, matcode);
		mm_write_mtx_crd_size(stdout, M, N, nz);
		for (i = 0; i < nz; i++)
			fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);

		double* result = (double *)malloc(M * N * sizeof(double));
		convertToNonCompressedMatrix(M, N, nz, I, J, val, result);
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

int main(int argc, char** argv) 
{
  
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
  int row, col, idx;
  
  generateVector(nrows, x);
  generateMatrix(nrows, ncols, A);

  // We check the matrix and the vector
  fprintf(stdout, "Matrix:\n");
  printMatrix(nrows, ncols, A);
  fprintf(stdout, "Vector:\n");
  printVector(nrows, x);
  
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
