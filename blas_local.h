#ifndef BLAS_LOCAL_H
#define BLAS_LOCAL_H

#include <complex>

namespace BLAS {

enum MATRIX_ORDER {
   RowMajor=101,
   ColMajor=102
};

enum MATRIX_TRANSPOSE {
   NoTrans=111,
   Trans=112,
   ConjTrans=113
};

enum MATRIX_UPLO {
	Upper=121,
	Lower=122
};

float nrm2(const unsigned int N, const float *X, const unsigned int incX);
double nrm2(const unsigned int N, const double *X, const unsigned int incX);
float nrm2(const unsigned int N, const std::complex<float> *X, const unsigned int incX);
double nrm2(const unsigned int N, const std::complex<double> *X, const unsigned int incX);

float dot(const unsigned int N, const float *X, const unsigned int incX, const float *Y, const unsigned int incY);
double dot(const unsigned int N, const double *X, const unsigned int incX, const double *Y, const unsigned int incY);

std::complex<float> dotu(const unsigned int N, const std::complex<float> *X, const unsigned int incX, const std::complex<float> *Y, const unsigned int incY);
std::complex<double> dotu(const unsigned int N, const std::complex<double> *X, const unsigned int incX, const std::complex<double> *Y, const unsigned int incY);

std::complex<float> dotc(const unsigned int N, const std::complex<float> *X, const unsigned int incX, const std::complex<float> *Y, const unsigned int incY);
std::complex<double> dotc(const unsigned int N, const std::complex<double> *X, const unsigned int incX, const std::complex<double> *Y, const unsigned int incY);

void copy(const int N, const float *X, const int incX, float *Y, const int incY);
void copy(const int N, const double *X, const int incX, double *Y, const int incY);
void copy(const int N, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY);
void copy(const int N, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY);

void scal(const int N, const float *alpha, float *X, const int incX);
void scal(const int N, const double *alpha, double *X, const int incX);
void scal(const int N, const std::complex<float> *alpha, std::complex<float> *X, const int incX);
void scal(const int N, const std::complex<double> *alpha, std::complex<double> *X, const int incX);

void axpy (const int N, const float *alpha, const float *X, const int incX, float *Y, const int incY);
void axpy (const int N, const double *alpha, const double *X, const int incX, double *Y, const int incY);
void axpy (const int N, const std::complex<float> *alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY);
void axpy (const int N, const std::complex<double> *alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY);

/*
void ger(const enum MATRIX_ORDER Order, const int M, const int N, const datatype alpha, const datatype *X, const int incX, const datatype *Y, const int incY, datatype *A, const int lda);
void geru(const enum MATRIX_ORDER Order, const int M, const int N, const cdatatype alpha, const cdatatype *X, const int incX, const cdatatype *Y, const int incY, cdatatype *A, const int lda);
void gerc(const enum MATRIX_ORDER Order, const int M,const int N, const cdatatype alpha, const cdatatype *X, const int incX, const cdatatype *Y, const int incY, cdatatype *A, const int lda);

*/

void her(const enum MATRIX_ORDER Order, const enum MATRIX_UPLO Uplo, const int N, const float alpha, const std::complex<float> *X, const int incX, std::complex<float> *A, const int lda);
void her(const enum MATRIX_ORDER Order, const enum MATRIX_UPLO Uplo, const int N, const double alpha, const std::complex<double> *X, const int incX, std::complex<double> *A, const int lda);

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *B, const int ldb, const std::complex<float> beta, std::complex<float> *C, const int ldc);

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *B, const int ldb, const std::complex<double> beta, std::complex<double> *C, const int ldc);


void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const float alpha, const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY);

void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY);

void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incX, const std::complex<float> beta, std::complex<float> *Y, const int incY);

void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incX, const std::complex<double> beta, std::complex<double> *Y, const int incY);


}

#endif
