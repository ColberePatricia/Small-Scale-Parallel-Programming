// This file contains the test functions

#ifndef TEST_H_
#define TEST_H_

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "matrixPreprocessing.h"
#include "matrixVector.h"


// Do all the unit tests for matrix processing
void testMatrixProcessing();

// These are the unit tests for matrix processing
void testIsEqual();
void testIsEqualInt();
void testBubbleSort();
void testMaximumMatrix();

// We do the tests of reorder on a matrix where it will modify I, J and val
void testReorderI();
void testReorderJ();
void testreorderVal();

void testConvertToNonCompressedMatrix();

void testGetCSR_IRP();
void testGetCSR_JA();
void testGetCSR_AS();

void testGetELLPACK_MAXNZ();
void testGetELLPACK_JA();
void testGetELLPACK_AS();

void testReadMatrix();



// Do all the unit tests for the matrix vector product
void testMatrixVectorProduct();

// These are the unit tests for the matrix vector product
void testMatrixVector();
void testMatrixVectorCSR();
void testMatrixVectorELLPACK();



// These functions will be used to test matrix processing and the product y<-Ax
/*
Matrix created
(11 12 0 0)
(0 22 23 0)
(0 0  33 0)
(0 0 43 44)
*/
int* createTestI();
int* createTestJ();
double* createTestVal();
/*
Matrix created
(1 7 0 0)
(0 2 8 0)
(5 0 3 9)
(0 6 0 4)
*/
int* createTestIReorder();
int* createTestJReorder();
double* createTestValReorder();
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
int* createTestIIRP();


#endif