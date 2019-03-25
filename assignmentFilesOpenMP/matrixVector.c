#include "matrixVector.h"


inline int min(int a, int b) { return a < b ? a : b; }

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
	int i, j;
	// We try to unroll to 4 if possible
#pragma omp parallel for shared(x,y,IRP, JA, AS) private(i,j,temp)
	for (int i = 0;i < M;i++) {
		temp = 0;
		if (IRP[i + 1] - IRP[i] >= 4) {
			for (int j = IRP[i];j <= (IRP[i + 1] - 1) - (IRP[i + 1] - 1) % 4;j += 4) {
				temp += AS[j] * x[JA[j]] + AS[j + 1] * x[JA[j + 1]] + AS[j + 2] * x[JA[j + 2]] + AS[j + 3] * x[JA[j + 3]];
			}
			for (int j = (IRP[i + 1] - 1) - (IRP[i + 1] - 1) % 4;j <= IRP[i + 1] - 1;j++) {
				temp += AS[j] * x[JA[j]];
			}
		}
		else {
			for (int j = IRP[i];j <= IRP[i + 1] - 1;j++) {
				temp += AS[j] * x[JA[j]];
			}
		}
		y[i] = temp;

	}
	return y;
}




// Implementation of matrix-vector product with the matrix in ELLPACK
double* MatrixVectorELLPACK(int M, int N, int MAXNZ, int* JA, double* AS, double* x, double* y) {
	double t, t0, t1, t2, t3;
	int i, j, idx; // idx is the index of (i,j)

	// We unroll to 4 to reduce the loading time
#pragma omp parallel for shared(x,y,MAXNZ, JA, AS,M,N) private(i,j,idx,t0, t1, t2) 
	for (i = 0;i < M - M % 4;i += 4) {
		t0 = 0;
		t1 = 0;
		t2 = 0;
		t3 = 0;

		for (j = 0;j < MAXNZ - MAXNZ % 2;j += 2) {
			t0 += AS[(i + 0)*MAXNZ + j + 0] * x[JA[(i + 0)*MAXNZ + j + 0]] + AS[(i + 0)*MAXNZ + j + 1] * x[JA[(i + 0)*MAXNZ + j + 1]];
			t1 += AS[(i + 1)*MAXNZ + j + 0] * x[JA[(i + 1)*MAXNZ + j + 0]] + AS[(i + 1)*MAXNZ + j + 1] * x[JA[(i + 1)*MAXNZ + j + 1]];
			t2 += AS[(i + 2)*MAXNZ + j + 0] * x[JA[(i + 2)*MAXNZ + j + 0]] + AS[(i + 2)*MAXNZ + j + 1] * x[JA[(i + 2)*MAXNZ + j + 1]];
			t3 += AS[(i + 3)*MAXNZ + j + 0] * x[JA[(i + 3)*MAXNZ + j + 0]] + AS[(i + 3)*MAXNZ + j + 1] * x[JA[(i + 3)*MAXNZ + j + 1]];
		}

		for (j = MAXNZ - MAXNZ % 2;j < MAXNZ;j++) {
			t0 += AS[(i + 0)*MAXNZ + j] * x[JA[(i + 0)*MAXNZ + j]];
			t1 += AS[(i + 1)*MAXNZ + j] * x[JA[(i + 1)*MAXNZ + j]];
			t2 += AS[(i + 2)*MAXNZ + j] * x[JA[(i + 2)*MAXNZ + j]];
			t3 += AS[(i + 3)*MAXNZ + j] * x[JA[(i + 3)*MAXNZ + j]];
		}
		y[i + 0] = t0;
		y[i + 1] = t1;
		y[i + 2] = t2;
		y[i + 3] = t3;
	}

	for (i = M - M % 4;i < M;i++) {
		t = 0.0;
		for (j = 0;j < MAXNZ;j++) {
			idx = i * MAXNZ + j;
			t += AS[idx] * x[JA[idx]];
		}
		y[i] = t;
	}

	return y;
}
