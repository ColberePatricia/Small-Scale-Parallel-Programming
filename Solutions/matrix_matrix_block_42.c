// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
// Multiplies  two matrices. Matrices are stored in linear memory in row-major order,
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"

// Matrix dimensions.
const int ROWS = 1024;
const int COLS = 1024;
const int ntimes = 5;
inline double dmin ( double a, double b ) { return a < b ? a : b; }
inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

#define BS  32
void innerMatrixMultiply(int nr,int nc,int nk,
			 double AA[BS][BS],
			 double BB[BS][BS],
			 double CC[BS][BS])
{
  for (int i=0; i<BS; i += 4) {
    for (int j=0; j<BS; j +=2 ) {
      double t00=CC[i+0][j+0];
      double t10=CC[i+1][j+0];
      double t20=CC[i+2][j+0];
      double t30=CC[i+3][j+0];
      double t01=CC[i+0][j+1];
      double t11=CC[i+1][j+1];
      double t21=CC[i+2][j+1];
      double t31=CC[i+3][j+1];

      for (int k=0; k<BS; k += 4) {
	double a00=AA[i+0][k+0];
	double a01=AA[i+0][k+1];
	double a02=AA[i+0][k+2];
	double a03=AA[i+0][k+3];
	double a10=AA[i+1][k+0];
	double a11=AA[i+1][k+1];
	double a12=AA[i+1][k+2];
	double a13=AA[i+1][k+3];
	double a20=AA[i+2][k+0];
	double a21=AA[i+2][k+1];
	double a22=AA[i+2][k+2];
	double a23=AA[i+2][k+3];
	double a30=AA[i+3][k+0];
	double a31=AA[i+3][k+1];
	double a32=AA[i+3][k+2];
	double a33=AA[i+3][k+3];
	double b00=BB[k+0][j+0];
	double b10=BB[k+1][j+0];
	double b20=BB[k+2][j+0];
	double b30=BB[k+3][j+0];
	double b01=BB[k+0][j+1];
	double b11=BB[k+1][j+1];
	double b21=BB[k+2][j+1];
	double b31=BB[k+3][j+1];


	t00 += a00*b00 + a01*b10 + a02*b20 + a03*b30;
	t10 += a10*b00 + a11*b10 + a12*b20 + a13*b30;
	t20 += a20*b00 + a21*b10 + a22*b20 + a23*b30;
	t30 += a30*b00 + a31*b10 + a32*b20 + a33*b30;

	t01 += a00*b01 + a01*b11 + a02*b21 + a03*b31;
	t11 += a10*b01 + a11*b11 + a12*b21 + a13*b31;
	t21 += a20*b01 + a21*b11 + a22*b21 + a23*b31;
	t31 += a30*b01 + a31*b11 + a32*b21 + a33*b31;

      }
    CC[i+0][j+0]=t00;
    CC[i+1][j+0]=t10;
    CC[i+2][j+0]=t20;
    CC[i+3][j+0]=t30;

    CC[i+0][j+1]=t01;
    CC[i+1][j+1]=t11;
    CC[i+2][j+1]=t21;
    CC[i+3][j+1]=t31;


    }
  }
}

void copyToBuffer(int nr, int nc, const double *X, int colsx, double XX[BS][BS])
{
  int i,j,idx;
  for ( i=0; i<nr; i++){
    for ( j=0; j<nc; j++){
      idx = (i) * colsx +j;
      XX[i][j] = X[idx];
    }
    for ( j=nc; j<BS; j++){
      XX[i][j] = 0.0;
    }
  }
  for ( i=0; i<nr; i++)
    for ( j=0; j<BS; j++)
      XX[i][j] = 0.0;
  
}

void copyFromBuffer(int nr, int nc, double *X, int colsx, double XX[BS][BS])
{
  int i,j,idx;
  for ( i=0; i<nr; i++){
    for ( j=0; j<nc; j++){
      idx = (i) * colsx +j;
      X[idx] = XX[i][j] ;
    }
  }
}

void zeroBuffer(double XX[BS][BS])
{
  int i,j;
  for ( i=0; i<BS; i++){
    for ( j=0; j<BS; j++){
      XX[i][j] =0.0;
    }
  }
}

// Block and unroll  implementation of matrix multiply
void MatrixMultiply(int rowsc, int colsc, int colsa,
		    const double* A, const double* B,
		    double* restrict C) 
{
  int rowc,colc,rowb,cola,idxa,idxb,idxc;
#pragma omp parallel for shared(A,B,C) private(rowc,colc,cola,idxa,idxb,idxc)
  for (rowc = 0; rowc < rowsc; rowc +=BS) {
    double AA[BS][BS],BB[BS][BS],CC[BS][BS];
    int nr=min(BS,rowsc-rowc);
    int nc, nk;
    for (colc = 0; colc < colsc; colc += BS) {
      nc=min(BS,colsc-colc);
      idxc = (rowc) * colsc + colc;
      //copyToBuffer(nr,nc,&(C[idxc]),colsc,CC);
      zeroBuffer(CC);      
      for (cola=0; cola < colsa; cola += BS) {
	nk=min(BS,colsa-cola);
	idxa = rowc * colsa + cola;
	copyToBuffer(nr,nk,&(A[idxa]),colsa,AA);
	idxb = cola * colsc + colc;
	copyToBuffer(nk,nc,&(B[idxb]),colsc,BB);
	innerMatrixMultiply(nr,nc,nk,AA,BB,CC);
      }
      copyFromBuffer(nr,nc,&(C[idxc]),colsc,CC);      
    }
  }
}

int main(int argc, char** argv) 
{
  if (argc < 4) {
    fprintf(stderr,"Usage: %s  rowsc colsc colsa\n",argv[0]);
  }
  int nrowsc=atoi(argv[1]);
  int ncolsc=atoi(argv[2]);
  int ncolsa=atoi(argv[3]);
  
  
  // ----------------------- Host memory initialisation ----------------------- //
  
  double* A = (double*) malloc(sizeof(double)*nrowsc*ncolsa);
  double* B = (double*) malloc(sizeof(double)*ncolsa*ncolsc);
  double* C = (double*) malloc(sizeof(double)*nrowsc*ncolsc);
  int row, col, idx;
  
  srand(12345);
  for ( row = 0; row < nrowsc; ++row) {
    for ( col = 0; col < ncolsa; ++col) {
      idx = row * ncolsa + col;
      A[idx] = 100.0f * ((double) rand()) / RAND_MAX;
    }
  }
  for ( row = 0; row < ncolsa; ++row) {
    for ( col = 0; col < ncolsc; ++col) {
      idx = row * ncolsc + col;
      B[idx] = 100.0f * ((double) rand()) / RAND_MAX;
    }
  }
  
  double tmlt = 1e100;
  for (int try=0; try < ntimes; try ++ ) {
    double t1 = wtime();
    MatrixMultiply(nrowsc,ncolsc,ncolsa, A, B, C);
    double t2 = wtime();
    tmlt = dmin(tmlt,(t2-t1));
  }
  double mflops = (2.0e-6)*nrowsc*ncolsc*ncolsa/tmlt;
#pragma omp parallel 
  {
#pragma omp master
    {
      fprintf(stdout,"Multiplying matrices of size %d x %d (%d) (block_42) with %d threads: time %lf  MFLOPS %lf \n",
	      nrowsc,ncolsc,ncolsa,omp_get_num_threads(),tmlt,mflops);
    }
  }
  free(A);
  free(B);
  free(C);
  return 0;
}
