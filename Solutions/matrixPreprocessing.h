// This file contains the functions to preprocess a matrix before we can apply our calculations on it

#ifndef MATRIXPREPROCESSING_H_
#define MATRIXPREPROCESSING_H_

#include <stdlib.h>
#include <stdio.h>
#include "mmio.h"
#include "print.h"


double* convertToNonCompressedMatrix(int M, int N, int nz, int *I, int *J, double* val);

// Returns the IRP of the CSR
double* getCSR_IRP(int M, int N, int nz, int *I, int *J, double* val);
// Returns the JA of the CSR
int* getCSR_JA(int M, int N, int nz, int *I, int *J, double* val);
// Returns the AS of the CSR
double* getCSR_AS(int M, int N, int nz, int *I, int *J, double* val);
// Returns the MAXNZ of the ELLPACK
int getELLPACK_MAXNZ(int M, int N, int nz, int *I, int *J, double* val);
// Returns the JA of the ELLPACK
int* getELLPACK_JA(int M, int N, int nz, int *I, int *J, double* val, int MAXNZ);
// Returns the AS of the ELLPACK
double* getELLPACK_AS(int M, int N, int nz, int *I, int *J, double* val, int MAXNZ);

// Read a matrix from the input
void readMatrix(char* fileName);


#endif

