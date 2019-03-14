// This file contains the functions to print a matrix or a vector

#ifndef PRINT_H_
#define PRINT_H_

#include <stdlib.h>
#include <stdio.h>

// Print the matrix in the console
void printMatrix(int rows, int cols, const double* A);

// Print the vector in the console
void printVector(int size, const double* vector);



#endif