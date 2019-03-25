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
// using one thread per row
__global__ void MatrixVectorCSR(int M, int* IRP, int* JA, double* AS, double* x, double* y) {
	double temp;
	int tr = threadIdx.x;
	int i = blockIdx.x*blockDim.x + tr;
	if (i < M) {
		temp = 0;
		for (int j = IRP[i];j <= IRP[i + 1] - 1;j++) {
			temp += AS[j] * x[JA[j]];
		}
		y[i] = temp;
	}
}




// Implementation of matrix-vector product with the matrix in ELLPACK
// using a block of threads for each block of rows.
__global__ void MatrixVectorELLPACK(int M, int N, int MAXNZ, int* JA, double* AS, double* x, double* y) {
	__shared__ double ax[16][64];
	double temp;
	int tr     = threadIdx.y;
  int tc     = threadIdx.x;
  int i    = blockIdx.x*blockDim.y + tr;
  ax[tr][tc] = 0.0;
	if (i < M) {
		int idx = i * MAXNZ + tc;
		temp = 0;
		int j;
		for (j = tc;j < MAXNZ;j+=64) {
			temp += AS[idx] * x[JA[idx]];
			idx+=64; // The size of JA is M * MAXNZ
		}
		if (j<MAXNZ){
			temp += AS[idx] * x[JA[idx]];
		}
		ax[tr][tc] = temp;
	}
	__syncthreads();
	for (int s=64/2; s >32; s >>=1){
		if (tc<s)
			ax[tr][tc] += ax[tr][tc+s];
			__syncthreads();
	}
	for (int s=min(32,64/2); s >0; s >>=1){
		if (tc<s)
			ax[tr][tc] += ax[tr][tc+s];
	}

	if ((tc == 0)&&(i<M))
		y[i] = ax[tr][tc];
}

