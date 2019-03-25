// This file contains the functions to preprocess a matrix before we can apply our calculations on it

#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers


// Implementation of matrix-vector product with the matrix not compressed
double* MatrixVector(int rows, int cols, const double* A, const double* x, double* y);
// Implementation of matrix-vector product with the matrix in CSR
// using a block of threads for each block of rows.
__global__ void MatrixVectorCSR(int M, int* IRP, int* JA, double* AS, double* x, double* y);
// Implementation of matrix-vector product with the matrix in ELLPACK
// using a block of threads for each block of rows. 
__global__ void MatrixVectorELLPACK(int M, int N, int MAXNZ, int* JA, double* AS, double* x, double* y);


#endif
