#include "matrixPreprocessing.h"


double* convertToNonCompressedMatrix(int M, int N, int nz, int *I, int *J, double* val, double* result) {
	// We reserve memory for the result matrix
	result = (double *)malloc(M * N * sizeof(double));
	// We initialize the matrix with zeroes
	for (int i = 0;i < M*N;i++)
		result[i] = 0;

	int idx;

	for (int i = 0;i < nz;i++) {
		idx = I[i] * N + J[i];
		//fprintf(stdout, "I=%d, J=%d, idx=%d\n", I[i], J[i], idx);
		result[idx] = val[i];
	}
	
	return result;
}

double* convertToELLPACK(int M, int N, int nz, int *I, int *J, double* val, double* result) {

}

double* convertToCSR(int M, int N, int nz, int *I, int *J, double* val, double* result) {

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

		
		// TO DO REMOVE!!!
		double* result = (double *)malloc(M * N * sizeof(double));
		result = convertToNonCompressedMatrix(M, N, nz, I, J, val, result);
		fprintf(stdout, "RESULT:\n\n");
		printMatrix(M, N, result);

	}
}