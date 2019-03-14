// This file contains the functions to preprocess a matrix before we can apply our calculations on it

#ifndef MATRIXPREPROCESSING_H_
#define MATRIXPREPROCESSING_H_

#include <stdlib.h>
#include <stdio.h>
#include "mmio.h"
#include "print.h"


double* convertToNonCompressedMatrix(int M, int N, int nz, int *I, int *J, double* val, double* result);
double* convertToELLPACK(int M, int N, int nz, int *I, int *J, double* val, double* result);
double* convertToCSR(int M, int N, int nz, int *I, int *J, double* val, double* result);

// Read a matrix from the input
void readMatrix(char* fileName);


#endif

