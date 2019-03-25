#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"
#include "mmio.h"
#include "matrixPreprocessing.h"
#include "test.h"
#include "matrixVector.h"


inline double dmin(double a, double b) { return a < b ? a : b; }


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
#pragma omp parallel
		{
#pragma omp master
			{
				fprintf(stderr, "Usage: %s [martix-market-filename]\nWe will be using the cage4 matrix as default\n", argv[0]);
				//fileName = "D:\\Cranfield work\\Small Scale Parallel Programming\\matrices\\cage4.mtx";
				fileName = "../matrices/cage4.mtx";
			}
		}
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


		double tmltCSR = 1e100;
		double tmltELLPACK = 1e100;
		// We do the product with CSR
		double t1CSR = wtime();
		for (int i=0;i<10;i++)
			MatrixVectorCSR(M, IRP, JA_CSR, AS_CSR, x, y);
		double t2CSR = wtime();
		tmltCSR = dmin(tmltCSR, (t2CSR - t1CSR));
		double mflopsCSR = (2.0e-6)*M*N*10 / tmltCSR;

		// We do the product for ELLPACK
		double t1ELLPACK = wtime();
		for (int i=0;i<10;i++)
			MatrixVectorELLPACK(M, N, MAXNZ, JA_ELLPACK, AS_ELLPACK, x, y);
		double t2ELLPACK = wtime();
		tmltELLPACK = dmin(tmltELLPACK, (t2ELLPACK - t1ELLPACK));
		double mflopsELLPACK = (2.0e-6)*M*N*10 / tmltELLPACK;


		// We print our results
#pragma omp parallel
		{
#pragma omp master
			{
				fprintf(stdout, "CSR: Matrix-Vector product of size %d x %d with %d threads: time %lf  MFLOPS %lf \n", M, N, omp_get_num_threads(), (tmltCSR/10), mflopsCSR);
				fprintf(stdout, "ELLPACK: Matrix-Vector product of size %d x %d with %d threads: time %lf  MFLOPS %lf \n", M, N, omp_get_num_threads(), (tmltELLPACK/10), mflopsELLPACK);
			}
		}

		// We free the matrices and vectors
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
