#include "print.h"

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