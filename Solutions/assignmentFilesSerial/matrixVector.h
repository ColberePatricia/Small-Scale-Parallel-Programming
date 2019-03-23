// This file contains the functions to preprocess a matrix before we can apply our calculations on it

#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_

#include <stdlib.h>
#include <stdio.h>


// Implementation of matrix-vector product with the matrix not compressed
double* MatrixVector(int rows, int cols, const double* A, const double* x, double* y);
// Implementation of matrix-vector product with the matrix in CSR
double* MatrixVectorCSR(int M, int* IRP, int* JA, double* AS, double* x, double* y);
// Implementation of matrix-vector product with the matrix in ELLPACK
double* MatrixVectorELLPACK(int M, int N, int MAXNZ, int* JA, double* AS, double* x, double* y);


#endif

