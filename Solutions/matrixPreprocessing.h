// This file contains the functions to preprocess a matrix before we can apply our calculations on it

#ifndef MATRIXPREPROCESSING_H_
#define MATRIXPREPROCESSING_H_

#include <stdlib.h>
#include <stdio.h>
#include "mmio.h"
#include "print.h"

// Return 0 if the matrices are different, 1 if they are the same
int isEqual(int size, double* matrix1, double* matrix2);
int isEqualInt(int size, int * matrix1, int * matrix2);
double* bubbleSort(int size, double* matrix);
// Returns the maximum for a matrix of int
int maximumMatrix(int size, int* matrix);

// reorder I, J and val so that the I is the increasing matrix instead of the J because we are coding in C
int* reorderI(int nz, int* I);
int* reorderJ(int nz, int* I, int* J);
double* reorderVal(int nz, int* I, double* val);

double* convertToNonCompressedMatrix(int M, int N, int nz, int *I, int *J, double* val);

// Returns the IRP of the CSR
int* getCSR_IRP(int M, int nz, int *I);
// Returns the JA of the CSR
int* getCSR_JA(int *J);
// Returns the AS of the CSR
double* getCSR_AS(double* val);

// Returns the MAXNZ of the ELLPACK
int getELLPACK_MAXNZ(int nz, int *I);
// Returns the JA of the ELLPACK
int* getELLPACK_JA(int M, int nz, int *I, int *J, int MAXNZ);
// Returns the AS of the ELLPACK
double* getELLPACK_AS(int M, int nz, int *I, double* val, int MAXNZ);

// Read a matrix from the input
void readMatrix(char* fileName);


#endif

