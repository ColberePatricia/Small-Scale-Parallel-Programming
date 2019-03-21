#include "matrixVector.h"

// Simple CPU implementation of matrix-vector product
double* MatrixVector(int rows, int cols, const double* A, const double* x, double* y) {
	int row, col, idx;
	double t;
	for (row = 0; row < rows; ++row) {
		t = 0.0;
		for (col = 0; col < cols; ++col) {
			idx = row * cols + col;
			t = t + A[idx] * x[col];
		}
		y[row] = t;
	}
	return y;
}

// Implementation of matrix-vector product with the matrix in CSR
double* MatrixVectorCSR(int M, int* IRP, int* JA, double* AS, double* x, double* y) {
	double temp;
	for (int i = 0;i < M;i++) {
		temp = 0;
		for (int j = IRP[i];j <= IRP[i + 1] - 1;j++) {
			//fprintf(stdout, "i: %d, j: %d\n", i, j);
			//fprintf(stdout, "IRP: %d\n", IRP[i]);
			temp += AS[j] * x[JA[j]];
		}
		y[i] = temp;
	}
	return y;
}

// Implementation of matrix-vector product with the matrix in ELLPACK
double* MatrixVectorELLPACK(int M, int N, int MAXNZ, int* JA, double* AS, double* x, double* y) {
	double temp;
	int idx; // The index of (i,j)
	for (int i = 0;i < M;i++) {
		temp = 0;
		for (int j = 0;j < MAXNZ;j++) {
			idx = i * MAXNZ + j; // The size of JA is M * MAXNZ
			temp += AS[idx] * x[JA[idx]];
		}
		y[i] = temp;
	}
	return y;
}