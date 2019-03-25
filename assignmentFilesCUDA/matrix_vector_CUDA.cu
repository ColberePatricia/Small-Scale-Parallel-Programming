#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers
#include "wtime.h"
#include "mmio.h"
#include "matrixPreprocessing.h"
#include "test.h"
#include "matrixVector.h"

//Simple dimension: define a 1D block structure
#define BD 256
const dim3 BLOCK_DIM(BD);

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




int main(int argc, char** argv)
{

	testMatrixProcessing();
	testMatrixVectorProduct();

	char* fileName;
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\nWe will be using the cage4 matrix as default\n", argv[0]);
		//fileName = "D:\\Cranfield work\\Small Scale Parallel Programming\\matrices\\cage4.mtx";
		fileName = "../matrices/cage4.mtx";
	}
	else {
		fileName = argv[1];
	}


	// We read the file of the matrix
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
		for (i = 0; i < nz; i++) {
			fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
			I[i]--;  // adjust from 1-based to 0-based
			J[i]--;
		}
		if (f != stdin) fclose(f);


		// We now have the matrix with I, J and val
	// We will generate our CSR matrix from them
		int* IRP = (int *)malloc((M + 1) * sizeof(int));
		int* JA_CSR = (int *)malloc(nz * sizeof(int));
		double* AS_CSR = (double *)malloc(nz * sizeof(double));
		IRP = getCSR_IRP(M, nz, I);
		JA_CSR = getCSR_JA(nz, I, J);
		AS_CSR = getCSR_AS(nz, I, val);

		// We will generate our ELLPACK matrix
		int MAXNZ;
		MAXNZ = getELLPACK_MAXNZ(nz, I);
		int* JA_ELLPACK = (int *)malloc(M * MAXNZ * sizeof(int));
		double* AS_ELLPACK = (double *)malloc(M * MAXNZ * sizeof(double));
		JA_ELLPACK = getELLPACK_JA(M, nz, I, J, MAXNZ);
		AS_ELLPACK = getELLPACK_AS(M, nz, I, val, MAXNZ);


		// We initiate our matrices for the product
		double* x = (double*)malloc(sizeof(double)*N);
		double* y = (double*)malloc(sizeof(double)*M);

		// We generate randomly x of size N
		generateVector(N, x);



		// We create our CUDA matrices
		double *d_AS_CSR, *d_AS_ELLPACK, *d_x, *d_y;
		int *d_IRP, *d_JA_CSR, *d_JA_ELLPACK;
		checkCudaErrors(cudaMalloc((void**)&d_AS_CSR, nz * sizeof(double)));
		checkCudaErrors(cudaMalloc((void**)&d_AS_ELLPACK, M * MAXNZ * sizeof(double)));
		checkCudaErrors(cudaMalloc((void**)&d_x, N * sizeof(double)));
		checkCudaErrors(cudaMalloc((void**)&d_y, M * sizeof(double)));
		checkCudaErrors(cudaMalloc((void**)&d_IRP, (M + 1) * sizeof(int)));
		checkCudaErrors(cudaMalloc((void**)&d_JA_CSR, nz * sizeof(int)));
		checkCudaErrors(cudaMalloc((void**)&d_JA_ELLPACK, M * MAXNZ * sizeof(int)));

		// Copy matrices from the host (CPU) to the device (GPU).
		checkCudaErrors(cudaMemcpy(d_AS_CSR, AS_CSR, nz * sizeof(double), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_AS_ELLPACK, AS_ELLPACK, M * MAXNZ * sizeof(double), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_x, x, N * sizeof(double), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_IRP, IRP, (M + 1) * sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_JA_CSR, JA_CSR, nz * sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_JA_ELLPACK, JA_ELLPACK, M * MAXNZ * sizeof(int), cudaMemcpyHostToDevice));


		// Calculate the dimension of the grid of blocks (1D) necessary to cover all rows.
		const dim3 GRID_DIM((M - 1 + BLOCK_DIM.x) / BLOCK_DIM.x, 1);
		double flopcnt = 2.e-6*M*N;

		// Create the CUDA SDK timer.
		StopWatchInterface* timer = 0;
		sdkCreateTimer(&timer);


		// We do the product with CSR
		timer->start();
		for (int i=0;i<10;i++)
			MatrixVectorCSR << <GRID_DIM, BLOCK_DIM >> > (M, d_IRP, d_JA_CSR, d_AS_CSR, d_x, d_y);
		checkCudaErrors(cudaDeviceSynchronize());
		timer->stop();
		double gpuflops = 10 * flopcnt / timer->getTime();

		// We print our results
		fprintf(stdout, "CSR: Matrix-Vector product of size %d x %d: time %lf  GFLOPS %lf \n", M, N, (timer->getTime())/10, gpuflops);


		// We do the product for ELLPACK
		timer->reset();
		timer->start();
		for (int i = 0;i < 10;i++)
			MatrixVectorELLPACK << <GRID_DIM, BLOCK_DIM >> > (M, N, MAXNZ, d_JA_ELLPACK, d_AS_ELLPACK, d_x, d_y);
		checkCudaErrors(cudaDeviceSynchronize());
		timer->stop();
		gpuflops = 10 * flopcnt / timer->getTime();

		// We print our results
		fprintf(stdout, "ELLPACK: Matrix-Vector product of size %d x %d: time %lf  GFLOPS %lf \n", M, N, (timer->getTime())/10, gpuflops);



		// We free the matrices and vectors

		delete timer;

		checkCudaErrors(cudaFree(d_IRP));
		checkCudaErrors(cudaFree(d_x));
		checkCudaErrors(cudaFree(d_JA_CSR));
		checkCudaErrors(cudaFree(d_AS_CSR));
		checkCudaErrors(cudaFree(d_JA_ELLPACK));
		checkCudaErrors(cudaFree(d_AS_ELLPACK));

		free(IRP);
		free(JA_CSR);
		free(AS_CSR);
		free(JA_ELLPACK);
		free(AS_ELLPACK);
		free(x);
		free(y);
	}

  return 0;
}
