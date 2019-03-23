#include "test.h"

// Do all the unit tests for matrix processing
void testMatrixProcessing() {
	fprintf(stdout, "We are starting the tests for matrix processing\n");

	testIsEqual();
	testIsEqualInt();
	testBubbleSort();
	testMaximumMatrix();

	testReorderI();
	testReorderJ();
	testreorderVal();

	testConvertToNonCompressedMatrix();

	testGetCSR_IRP();
	testGetCSR_JA();
	testGetCSR_AS();

	testGetELLPACK_MAXNZ();
	testGetELLPACK_JA();
	testGetELLPACK_AS();

	testReadMatrix();

	fprintf(stdout, "The tests for matrix processing have all been successful!\n");
}


// These are the unit tests for matrix processing

void testIsEqual() {
	double* m1 = (double *)malloc(5 * sizeof(double));
	double* m2 = (double *)malloc(5 * sizeof(double));
	double* m3 = (double *)malloc(5 * sizeof(double));

	for (int i = 0;i < 5;i++) {
		m1[i] = m3[i] = i;
		m2[i] = i * 2 + 1;
	}

	// These matrices should be different
	int resDiff = isEqual(5, m1, m2);
	// These matrices should be equal
	int resEq = isEqual(5, m1, m3);

	// We check our results
	assert(resDiff == 0);
	assert(resEq == 1);
}

void testIsEqualInt() {
	int* m1 = (int *)malloc(5 * sizeof(int));
	int* m2 = (int *)malloc(5 * sizeof(int));
	int* m3 = (int *)malloc(5 * sizeof(int));

	for (int i = 0;i < 5;i++) {
		m1[i] = m3[i] = i;
		m2[i] = i * 2 + 1;
	}

	// These matrices should be different
	int resDiff = isEqualInt(5, m1, m2);
	// These matrices should be equal
	int resEq = isEqualInt(5, m1, m3);

	// We check our results
	assert(resDiff == 0);
	assert(resEq == 1);
}

void testBubbleSort() {
	double* toBeSorted = (double *)malloc(5 * sizeof(double));
	double* sorted = (double *)malloc(5 * sizeof(double));

	// We initialise the two arrays so that sorted is the result we expect for a toBeSorted once it will have been sorted
	for (int i = 0;i < 5;i++) {
		sorted[i] = i;
	}
	toBeSorted[0] = 2;
	toBeSorted[1] = 0;
	toBeSorted[2] = 4;
	toBeSorted[3] = 1;
	toBeSorted[4] = 3;


	toBeSorted = bubbleSort(5, toBeSorted);
	
	assert(isEqual(5, toBeSorted, sorted) == 1);
}
void testMaximumMatrix() {
	int* m1 = (int *)malloc(5 * sizeof(int));

	m1[0] = 2;
	m1[1] = 0;
	m1[2] = -4;
	m1[3] = -1;
	m1[4] = 3;

	assert(maximumMatrix(5, m1) == 3);
}

// We do the tests of reorder on a matrix where it will modify I, J and val
void testReorderI() {
	int nz = 9;
	int* I = (int *)malloc(nz * sizeof(int));
	int* resultShouldBe = (int *)malloc(nz * sizeof(int));
	resultShouldBe[0] = resultShouldBe[1] = 0;
	resultShouldBe[2] = resultShouldBe[3] = 1;
	resultShouldBe[4] = resultShouldBe[5] = resultShouldBe[6] = 2;
	resultShouldBe[7] = resultShouldBe[8] = 3;

	I = createTestIReorder();
	I = reorderI(nz, I);
	//printMatrixInt(1, nz, I);
	//printMatrixInt(1, nz, resultShouldBe);

	assert(isEqualInt(nz, I, resultShouldBe) == 1);

}
void testReorderJ() {
	int nz = 9;
	int* I = (int *)malloc(nz * sizeof(int));
	int* J = (int *)malloc(nz * sizeof(int));
	int* resultShouldBe = (int *)malloc(nz * sizeof(int));
	resultShouldBe[0] = resultShouldBe[4] = 0;
	resultShouldBe[1] = resultShouldBe[2] = resultShouldBe[7] = 1;
	resultShouldBe[3] = resultShouldBe[5] = 2;
	resultShouldBe[6] = resultShouldBe[8] = 3;

	I = createTestIReorder();
	J = createTestJReorder();
	J = reorderJ(nz, I, J);
	//printMatrixInt(1, nz, J);
	//printMatrixInt(1, nz, resultShouldBe);

	assert(isEqualInt(nz, J, resultShouldBe) == 1);

}
void testreorderVal() {
	int nz = 9;
	int* I = (int *)malloc(nz * sizeof(int));
	double* val = (double *)malloc(nz * sizeof(double));
	double* resultShouldBe = (double *)malloc(nz * sizeof(double));
	resultShouldBe[0] = 1;
	resultShouldBe[1] = 7;
	resultShouldBe[2] = 2;
	resultShouldBe[3] = 8;
	resultShouldBe[4] = 5;
	resultShouldBe[5] = 3;
	resultShouldBe[6] = 9;
	resultShouldBe[7] = 6;
	resultShouldBe[8] = 4;

	I = createTestIReorder();
	val = createTestValReorder();
	val = reorderVal(nz, I, val);
	//printMatrix(1, nz, val);
	//printMatrix(1, nz, resultShouldBe);

	assert(isEqual(nz, val, resultShouldBe) == 1);
}

void testConvertToNonCompressedMatrix() {
	// We check that the test matrix in the assignment pdf is converted as expected

	// We initialize every value obtained from reading the matrix market file
	int M, N, nz;
	M = N = 4;
	nz = 7;
	int* I = (int *)malloc(nz * sizeof(int));
	int* J = (int *)malloc(nz * sizeof(int));
	double* val = (double *)malloc(nz * sizeof(double));
	// result will contain the matrix once it is uncompressed to ckeck that we read it properly
	double* result = (double *)malloc(M * N * sizeof(double));
	double* resultShouldBe = (double *)malloc(M * N * sizeof(double));

	I = createTestI();
	J = createTestJ();
	val = createTestVal();

	for (int i = 0;i < M*N;i++)
		resultShouldBe[i] = 0;
	resultShouldBe[0] = 11;
	resultShouldBe[1] = 12;
	resultShouldBe[5] = 22;
	resultShouldBe[6] = 23;
	resultShouldBe[10] = 33;
	resultShouldBe[14] = 43;
	resultShouldBe[15] = 44;

	result = convertToNonCompressedMatrix(M, N, nz, I, J, val);

	assert(isEqual(M*N, result, resultShouldBe) == 1);
}

void testGetCSR_IRP() {
	// We test IRP on the matrix
	/*
	(0 0 0 0)
	(1 2 8 0)
	(5 0 3 0)
	(0 0 0 0)
	(0 0 0 0)
	(0 0 0 7)
	(0 6 0 0)
	*/
	// IRP should be [0, 0, 3, 5, 5, 5, 6, 7]

	int M = 7;
	int nz = 7;
	int* I = (int *)malloc(nz * sizeof(int));
	int* IRPCalc = (int *)malloc((M + 1) * sizeof(int));
	int* resultShouldBe = (int *)malloc((M + 1) * sizeof(int));
	I = createTestIIRP();

	resultShouldBe[0] = 0;
	resultShouldBe[1] = 0;
	resultShouldBe[2] = 3;
	resultShouldBe[3] = 5;
	resultShouldBe[4] = 5;
	resultShouldBe[5] = 5;
	resultShouldBe[6] = 6;
	resultShouldBe[7] = 7;

	IRPCalc = getCSR_IRP(M, nz, I);

	//printMatrixInt(1, M+1, IRPCalc);
	//printMatrixInt(1, M+1, resultShouldBe);
	assert(isEqualInt(M + 1, IRPCalc, resultShouldBe));

}
void testGetCSR_JA() {
	// We test JA on the example on the assignment pdf
	int nz = 7;
	int* I = (int *)malloc(nz * sizeof(int));
	int* JA = (int *)malloc(nz * sizeof(int));
	int* resultShouldBe = (int *)malloc(nz * sizeof(int));

	I = createTestI();
	JA = createTestJ();
	resultShouldBe[0] = 0;
	resultShouldBe[1] = 1;
	resultShouldBe[2] = 1;
	resultShouldBe[3] = 2;
	resultShouldBe[4] = 2;
	resultShouldBe[5] = 2;
	resultShouldBe[6] = 3;

	JA = getCSR_JA(nz, I, JA);

	//printMatrixInt(1, nz, J);
	//printMatrixInt(1, nz, resultShouldBe);
	assert(isEqualInt(nz, JA, resultShouldBe) == 1);


}
void testGetCSR_AS() {
	// We test AS on the example on the assignment pdf
	int nz = 7;
	double* AS = (double *)malloc(nz * sizeof(double));
	int* I = (int *)malloc(nz * sizeof(int));
	double* resultShouldBe = (double *)malloc(nz * sizeof(double));

	AS = createTestVal();
	I = createTestI();
	resultShouldBe[0] = 11;
	resultShouldBe[1] = 12;
	resultShouldBe[2] = 22;
	resultShouldBe[3] = 23;
	resultShouldBe[4] = 33;
	resultShouldBe[5] = 43;
	resultShouldBe[6] = 44;

	AS = getCSR_AS(nz, I, AS);

	//printMatrix(1, nz, AS);
	//printMatrix(1, nz, resultShouldBe);
	assert(isEqual(nz, AS, resultShouldBe) == 1);
}

void testGetELLPACK_MAXNZ() {
	// We test MAXNZ on the example on the assignment pdf
	int nz = 7;
	int* I = (int *)malloc(nz * sizeof(int));
	int MAXNZ;
	int resultShouldBe = 2;

	I = createTestI();

	MAXNZ = getELLPACK_MAXNZ(nz, I);
	//fprintf(stdout, "MAXNZ: %d\n", MAXNZ);
	assert(MAXNZ == resultShouldBe);
}
void testGetELLPACK_JA() {
	int nz = 7;
	int M = 4;
	int MAXNZ = 2;
	int* I = (int *)malloc(nz * sizeof(int));
	int* J = (int *)malloc(nz * sizeof(int));
	int* JA = (int *)malloc(M * MAXNZ * sizeof(int));
	int* resultShouldBe = (int *)malloc(M * MAXNZ * sizeof(int));

	I = createTestI();
	J = createTestJ();

	JA = getELLPACK_JA(M, nz, I, J, MAXNZ);
	resultShouldBe[0] = 0;
	resultShouldBe[1] = 1;
	resultShouldBe[2] = 1;
	resultShouldBe[3] = 2;
	resultShouldBe[4] = 2;
	resultShouldBe[5] = 2;
	resultShouldBe[6] = 2;
	resultShouldBe[7] = 3;

	//printMatrixInt(M, MAXNZ, JA);
	//printMatrixInt(M, MAXNZ, resultShouldBe);
	assert(isEqualInt(M*MAXNZ, JA, resultShouldBe));
}
void testGetELLPACK_AS() {
	int nz = 7;
	int M = 4;
	int MAXNZ = 2;
	int* I = (int *)malloc(nz * sizeof(int));
	double* val = (double *)malloc(nz * sizeof(double));
	double* AS = (double *)malloc(M * MAXNZ * sizeof(double));
	double* resultShouldBe = (double *)malloc(M * MAXNZ * sizeof(double));

	I = createTestI();
	val = createTestVal();

	AS = getELLPACK_AS(M, nz, I, val, MAXNZ);
	resultShouldBe[0] = 11;
	resultShouldBe[1] = 12;
	resultShouldBe[2] = 22;
	resultShouldBe[3] = 23;
	resultShouldBe[4] = 33;
	resultShouldBe[5] = 0;
	resultShouldBe[6] = 43;
	resultShouldBe[7] = 44;

	//printMatrix(M, MAXNZ, AS);
	//printMatrix(M, MAXNZ, resultShouldBe);
	assert(isEqual(M*MAXNZ, AS, resultShouldBe));
}

void testReadMatrix() {
	//char* matrixFile = "D:\\Cranfield work\\Small Scale Parallel Programming\\matrices\\cage4.mtx";
	char* matrixFile = "../matrices/cage4.mtx";
	readMatrix(matrixFile);
}



// Do all the unit tests for the matrix vector product
void testMatrixVectorProduct() {
	fprintf(stdout, "We are starting the tests for the matrix vector product\n");
	testMatrixVector();
	testMatrixVectorCSR();
	testMatrixVectorELLPACK();
	fprintf(stdout, "The tests for the matrix vector product have all been successful!\n\n");
}

// These are the unit tests for the matrix vector product
void testMatrixVector() {
	int rows, cols;
	rows = cols = 4;
	double* A = (double*)malloc(sizeof(double)*rows * cols);
	double* x = (double*)malloc(sizeof(double)*rows);
	double* y = (double*)malloc(sizeof(double)*rows);
	double* resultShouldBe = (double*)malloc(sizeof(double)*rows);

	for (int i = 0;i < rows*cols;i++)
		A[i] = 0;
	A[0] = 11;
	A[1] = 12;
	A[5] = 22;
	A[6] = 23;
	A[10] = 33;
	A[14] = 43;
	A[15] = 44;

	for (int i = 0;i < rows;i++)
		x[i] = i + 1;

	resultShouldBe[0] = 35;
	resultShouldBe[1] = 113;
	resultShouldBe[2] = 99;
	resultShouldBe[3] = 305;

	y = MatrixVector(rows, cols, A, x, y);

	//printMatrix(1, rows, y);
	//printMatrix(1, rows, resultShouldBe);

	assert(isEqual(rows, y, resultShouldBe) == 1);
}

void testMatrixVectorCSR() {
	int M = 4;
	int nz = 7;
	int* IRP = (int *)malloc((M + 1) * sizeof(int));
	int* JA = (int *)malloc(nz * sizeof(int));
	double* AS = (double *)malloc(nz * sizeof(double));
	double* x = (double*)malloc(sizeof(double)*M);
	double* y = (double*)malloc(sizeof(double)*M);
	double* resultShouldBe = (double*)malloc(sizeof(double)*M);

	IRP[0] = 0;
	IRP[1] = 2;
	IRP[2] = 4;
	IRP[3] = 5;
	IRP[4] = 7;

	JA[0] = 0;
	JA[1] = JA[2] = 1;
	JA[3] = JA[4] = JA[5] = 2;
	JA[6] = 3;

	AS[0] = 11;
	AS[1] = 12;
	AS[2] = 22;
	AS[3] = 23;
	AS[4] = 33;
	AS[5] = 43;
	AS[6] = 44;

	for (int i = 0;i < M;i++)
		x[i] = i + 1;

	resultShouldBe[0] = 35;
	resultShouldBe[1] = 113;
	resultShouldBe[2] = 99;
	resultShouldBe[3] = 305;

	y = MatrixVectorCSR(M, IRP, JA, AS, x, y);

	//printMatrix(1, M, y);
	//printMatrix(1, M, resultShouldBe);

	assert(isEqual(M, y, resultShouldBe) == 1);



	// We do a second test with the matrix with which we tested IRP
	M = 7;
	int N = 4;
	nz = 7;
	int* IRP2 = (int *)malloc((M + 1) * sizeof(int));
	int* JA2 = (int *)malloc(nz * sizeof(int));
	double* AS2 = (double *)malloc(nz * sizeof(double));
	double* x2 = (double*)malloc(sizeof(double)*N);
	double* y2 = (double*)malloc(sizeof(double)*M);
	double* resultShouldBe2 = (double*)malloc(sizeof(double)*M);


	IRP2[0] = 0;
	IRP2[1] = 0;
	IRP2[2] = 3;
	IRP2[3] = 5;
	IRP2[4] = 5;
	IRP2[5] = 5;
	IRP2[6] = 6;
	IRP2[7] = 7;

	// JA is reordered because we are in C so we calculate by row
	JA2[0] = JA2[3] = 0;
	JA2[1] = JA2[6] = 1;
	JA2[2] = JA2[4] = 2;
	JA2[5] = 3;

	AS2[0] = 1;
	AS2[1] = 2;
	AS2[2] = 8;
	AS2[3] = 5;
	AS2[4] = 3;
	AS2[5] = 7;
	AS2[6] = 6;

	for (int i = 0;i < N;i++)
		x2[i] = i + 1;

	resultShouldBe2[0] = 0;
	resultShouldBe2[1] = 29;
	resultShouldBe2[2] = 14;
	resultShouldBe2[3] = 0;
	resultShouldBe2[4] = 0;
	resultShouldBe2[5] = 28;
	resultShouldBe2[6] = 12;

	y2 = MatrixVectorCSR(M, IRP2, JA2, AS2, x2, y2);

	//printMatrix(1, M, y2);
	//printMatrix(1, M, resultShouldBe2);

	assert(isEqual(M, y2, resultShouldBe2) == 1);
}

void testMatrixVectorELLPACK() {
	//int M, int N, int MAXNZ, int* JA, double* AS, double* x, double* y
	int M, N;
	M = N = 4;
	int MAXNZ = 2;
	int* JA = (int *)malloc(M * MAXNZ * sizeof(int));
	double* AS = (double *)malloc(M * MAXNZ * sizeof(double));
	double* x = (double*)malloc(sizeof(double)*M);
	double* y = (double*)malloc(sizeof(double)*M);
	double* resultShouldBe = (double*)malloc(sizeof(double)*M);

	JA[0] = 0;
	JA[1] = JA[2] = 1;
	JA[3] = JA[4] = JA[5] = JA[6] = 2;
	JA[7] = 3;

	AS[0] = 11;
	AS[1] = 12;
	AS[2] = 22;
	AS[3] = 23;
	AS[4] = 33;
	AS[5] = 0;
	AS[6] = 43;
	AS[7] = 44;

	for (int i = 0;i < M;i++)
		x[i] = i + 1;

	resultShouldBe[0] = 35;
	resultShouldBe[1] = 113;
	resultShouldBe[2] = 99;
	resultShouldBe[3] = 305;

	y = MatrixVectorELLPACK(M, N, MAXNZ, JA, AS, x, y);

	//printMatrix(1, M, y);
	//printMatrix(1, M, resultShouldBe);

	assert(isEqual(M, y, resultShouldBe) == 1);

}




/*
Matrix created
(11 12 0 0)
(0 22 23 0)
(0 0  33 0)
(0 0 43 44)
*/

int* createTestI() {
	int nz = 7;
	int* I = (int *)malloc(nz * sizeof(int));
	I[0] = I[1] = 0;
	I[2] = I[3] = 1;
	I[4] = 2;
	I[5] = I[6] = 3;
	return I;
}

int* createTestJ() {
	int nz = 7;
	int* J = (int *)malloc(nz * sizeof(int));
	J[0] = 0;
	J[1] = J[2] = 1;
	J[3] = J[4] = J[5] = 2;
	J[6] = 3;
	return J;
}

double* createTestVal() {
	int nz = 7;
	double* val = (double *)malloc(nz * sizeof(double));
	val[0] = 11;
	val[1] = 12;
	val[2] = 22;
	val[3] = 23;
	val[4] = 33;
	val[5] = 43;
	val[6] = 44;
	return val;
}


/*
Matrix created
(1 7 0 0)
(0 2 8 0)
(5 0 3 9)
(0 6 0 4)
*/

int* createTestIReorder() {
	int nz = 9;
	int* I = (int *)malloc(nz * sizeof(int));
	I[0] = I[2] = 0;
	I[3] = I[5] = 1;
	I[1] = I[6] = I[7] = 2;
	I[4] = I[8] = 3;
	return I;
}

int* createTestJReorder() {
	int nz = 9;
	int* J = (int *)malloc(nz * sizeof(int));
	J[0] = J[1] = 0;
	J[2] = J[3] = J[4] = 1;
	J[5] = J[6] = 2;
	J[7] = J[8] = 3;
	return J;
}

double* createTestValReorder() {
	int nz = 9;
	double* val = (double *)malloc(nz * sizeof(double));
	val[0] = 1;
	val[1] = 5;
	val[2] = 7;
	val[3] = 2;
	val[4] = 6;
	val[5] = 8;
	val[6] = 3;
	val[7] = 9;
	val[8] = 4;
	return val;
}

/*
Matrix created
(0 0 0 0)
(1 2 8 0)
(5 0 3 0)
(0 0 0 0)
(0 0 0 0)
(0 0 0 7)
(0 6 0 0)
*/
int* createTestIIRP() {
	int nz = 7;
	int* I = (int *)malloc(nz * sizeof(int));
	I[0] = I[2] = I[4] = 1;
	I[1] = I[5] = 2;
	I[6] = 5;
	I[3] = 6;
	return I;
}
