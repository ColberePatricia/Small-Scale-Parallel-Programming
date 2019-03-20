#include "matrixPreprocessing.h"

double* bubbleSort(int size, double* matrix) {
	int i, j;
	for (i = 0; i < size - 1; ++i) {
		for (j = 0; j < size - 1 - i; ++j) {
			if (matrix[j] > matrix[j + 1])
				matrix[j + 1], matrix[j] = matrix[j], matrix[j + 1];
		}
	}
}

int maximumMatrix(int size, double* matrix) {
	int maximum = matrix[0];

	for (int i = 1; i < size; i++) {
		if (matrix[i] > maximum)
			maximum = matrix[i];
	}
	return maximum;
}

double* convertToNonCompressedMatrix(int M, int N, int nz, int *I, int *J, double* val) {
	// We reserve memory for the result matrix
	double* result = (double *)malloc(M * N * sizeof(double));
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

// Returns the IRP of the CSR
double* getCSR_IRP(int M, int N, int nz, int *I, int *J, double* val) {
	int* IRP = (int *)malloc((M+1) * sizeof(int));
	int k = 1; // k will be the current index+1 of the IRP array
	int i; // i will be the current index+1 of the I array
	int irp = 1; // irp will correspond to the value of IRP

	// We put the first values to zero if the first rows are empty
	while (I[0] != k) {
		IRP[k - 1] = 0;
		k++;
	}

	// Now I[0] == k
	IRP[k - 1] = irp; // The first value is always one
	k++;
	irp++;

	// We go through the I array
	for (i = 2;i <= nz;i++) {
		// If we are on a new row, we can put a new value in IRP
		if (I[i - 1] == I[i - 2] + 1) {
			IRP[k - 1] = irp;
			k++;
			irp++;
		} // We have skipped at least a row
		else if (I[i - 1] > I[i - 2] + 1) {
			// We need to input the previous value again as many times as there are skipped rows
			for (int skipped = 1;skipped <= I[i - 1] - I[i - 2] - 1; skipped++) {
				IRP[k - 1] = IRP[k - 2];
				k++;
			}
			// We also need to input the value corresponding to the new row
			IRP[k - 1] = irp;
			irp++;
		}
		else {
			// The value increases because we have stayed on the same row but are moving in the index of non zero values
			irp++;
		}
	}

	// The last value is the number of non zero values plus one
	IRP[M] = nz + 1;

	return IRP;
}

// Returns the JA of the CSR
int* getCSR_JA(int M, int N, int nz, int *I, int *J, double* val) {
	return J;
}

// Returns the AS of the CSR
double* getCSR_AS(int M, int N, int nz, int *I, int *J, double* val) {
	return val;
}

// Returns the MAXNZ of the ELLPACK
int getELLPACK_MAXNZ(int M, int N, int nz, int *I, int *J, double* val) {
	// We create an array that will contain the number of non zero for each row
	// from this array we will get the max, that is MAXNZ
	int* temp = (int *)malloc(nz * sizeof(int));
	// We initialise its values to zero

	for (int i = 0;i < nz;i++)
		temp[I[i]]++;

	return maximumMatrix(nz, temp);
}

// Returns the JA of the ELLPACK
int* getELLPACK_JA(int M, int N, int nz, int *I, int *J, double* val, int MAXNZ) {
	int* JA = (int *)malloc(M * MAXNZ * sizeof(int));
	int k, p, q;
	k = 1;
	int idx;

	for (p = 1;p <= M;p++) {
		for (q = 1;q <= MAXNZ;q++) {
			idx = (p - 1)*MAXNZ + (q - 1);
			fprintf(stdout, "p-1=%d, q-1=%d, idx=%d\n", p - 1, q - 1, idx);
			if (I[k - 1] == p) {
				JA[idx] = J[k - 1];
				k++;
			}
			else
				JA[idx] = J[k - 2];
		}
	}
	return JA;
}

// Returns the AS of the ELLPACK
double* getELLPACK_AS(int M, int N, int nz, int *I, int *J, double* val, int MAXNZ) {
	double* AS = (double *)malloc(M * MAXNZ * sizeof(double));
	int k, p, q;
	k = 1;
	int idx;

	for (p = 1;p <= M;p++) {
		for (q = 1;q <= MAXNZ;q++) {
			idx = (p - 1)*MAXNZ + (q - 1);
			fprintf(stdout, "p-1=%d, q-1=%d, idx=%d\n", p - 1, q - 1, idx);
			if (I[k - 1] == p) {
				AS[idx] = val[k - 1];
				k++;
			} else
				AS[idx] = 0;
		}
	}
	return AS;
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
		result = convertToNonCompressedMatrix(M, N, nz, I, J, val);
		fprintf(stdout, "RESULT:\n\n");
		printMatrix(M, N, result);

	}
}