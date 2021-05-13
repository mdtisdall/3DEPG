#include "blas_local.h"

#include "mkl.h" 

#include <complex>

using namespace std;

namespace BLAS {

float nrm2(const unsigned int N, const float *X, const unsigned int incX) {
	return cblas_snrm2(N, X, incX);
}

double nrm2(const unsigned int N, const double *X, const unsigned int incX) {
	return cblas_dnrm2(N, X, incX);
}

float nrm2(const unsigned int N, const complex<float> *X, const unsigned int incX) {
	return cblas_scnrm2(N, X, incX);
}

double nrm2(const unsigned int N, const complex<double> *X, const unsigned int incX) {
	return cblas_dznrm2(N, X, incX);
}


void scal(const int N, const float *alpha, float *X, const int incX) {
	cblas_sscal(N, *alpha, X, incX);
}

void scal(const int N, const double *alpha, double *X, const int incX) {
	cblas_dscal(N, *alpha, X, incX);
}

void scal(const int N, const complex<float> *alpha, complex<float> *X, const int incX) {
	cblas_cscal(N, alpha, X, incX);
}

void scal(const int N, const complex<double> *alpha, complex<double> *X, const int incX) {
	cblas_zscal(N, alpha, X, incX);
}


float dot(const unsigned int N, const float *X, const unsigned int incX, const float *Y, const unsigned int incY) {
  MKL_INT innerN = N;
  MKL_INT innerIncX = incX;
  MKL_INT innerIncY = incY;

	return cblas_sdot(innerN, X, innerIncX, Y, innerIncY);
}

double dot(const unsigned int N, const double *X, const unsigned int incX, const double *Y, const unsigned int incY) {
  MKL_INT innerN = N;
  MKL_INT innerIncX = incX;
  MKL_INT innerIncY = incY;

	return cblas_ddot(innerN, X, innerIncX, Y, innerIncY);
}


complex<float> dotu(const unsigned int N, const complex<float> *X, const unsigned int incX, const complex<float> *Y, const unsigned int incY) {
	complex<float> ret;
	
	cblas_cdotu_sub(N, X, incX, Y, incY, &ret);

	return ret;
}

complex<double> dotu(const unsigned int N, const complex<double> *X, const unsigned int incX, const complex<double> *Y, const unsigned int incY) {
	complex<double> ret;
	
	cblas_zdotu_sub(N, X, incX, Y, incY, &ret);

	return ret;
}


complex<float> dotc(const unsigned int N, const complex<float> *X, const unsigned int incX, const complex<float> *Y, const unsigned int incY) {
	complex<float> ret;
	
	cblas_cdotc_sub(N, X, incX, Y, incY, &ret);

	return ret;
}

complex<double> dotc(const unsigned int N, const complex<double> *X, const unsigned int incX, const complex<double> *Y, const unsigned int incY) {
	complex<double> ret;
	
	cblas_zdotc_sub(N, X, incX, Y, incY, &ret);

	return ret;
}


void copy(const int N, const float *X, const int incX, float *Y, const int incY) {
	cblas_scopy(N, X, incX, Y, incY);
}

void copy(const int N, const double *X, const int incX, double *Y, const int incY) {
	cblas_dcopy(N, X, incX, Y, incY);
}

void copy(const int N, const complex<float> *X, const int incX, complex<float> *Y, const int incY) {
	cblas_ccopy(N, X, incX, Y, incY);
}

void copy(const int N, const complex<double> *X, const int incX, complex<double> *Y, const int incY) {
	cblas_zcopy(N, X, incX, Y, incY);
}


void axpy (const int N, const float *alpha, const float *X, const int incX, float *Y, const int incY) {
	cblas_saxpy(N, *alpha, X, incX, Y, incY);
}

void axpy (const int N, const double *alpha, const double *X, const int incX, double *Y, const int incY) {
	cblas_daxpy(N, *alpha, X, incX, Y, incY);
}

void axpy (const int N,
    const complex<float> *alpha, const complex<float> *X, const int incX,
    complex<float> *Y, const int incY) {
	cblas_caxpy(N,
    (MKL_Complex8*)alpha, (MKL_Complex8*)X, incX,
    (MKL_Complex8*)Y, incY);
}

void axpy (const int N,
    const complex<double> *alpha, const complex<double> *X, const int incX,
    complex<double> *Y, const int incY) {
	cblas_zaxpy(N, (MKL_Complex16*)alpha, (MKL_Complex16*)X, incX,
    (MKL_Complex16*)Y, incY);
}

/*
void ger(const enum MATRIX_ORDER Order, const int M, const int N, const datatype alpha, const datatype *X, const int incX, const datatype *Y, const int incY, datatype *A, const int lda) {
	#ifdef SINGLE_PRECISION
	cblas_sger ((CBLAS_ORDER) Order, M, N, alpha, X, incX, Y, incY, A, lda);
	#elif DOUBLE_PRECISION
	cblas_dger ((CBLAS_ORDER) Order, M, N, alpha, X, incX, Y, incY, A, lda);
	#endif		
}

void geru(const enum MATRIX_ORDER Order, const int M, const int N, const cdatatype alpha, const cdatatype *X, const int incX, const cdatatype *Y, const int incY, cdatatype *A, const int lda) {
	#ifdef SINGLE_PRECISION
	cblas_cgeru ((CBLAS_ORDER) Order, M, N, &alpha, X, incX, Y, incY, A, lda);
	#elif DOUBLE_PRECISION
	cblas_zgeru ((CBLAS_ORDER) Order, M, N, &alpha, X, incX, Y, incY, A, lda);
	#endif		
}

void gerc(const enum MATRIX_ORDER Order, const int M,const int N, const cdatatype alpha, const cdatatype *X, const int incX, const cdatatype *Y, const int incY, cdatatype *A, const int lda) {
	#ifdef SINGLE_PRECISION
	cblas_cgerc ((CBLAS_ORDER) Order, M, N, &alpha, X, incX, Y, incY, A, lda);
	#elif DOUBLE_PRECISION
	cblas_zgerc ((CBLAS_ORDER) Order, M, N, &alpha, X, incX, Y, incY, A, lda);
	#endif		
}


*/
void her(const enum MATRIX_ORDER Order, const enum MATRIX_UPLO Uplo, const int N, const float alpha, const std::complex<float> *X, const int incX, std::complex<float> *A, const int lda) {
	cblas_cher ((CBLAS_ORDER) Order, (CBLAS_UPLO) Uplo, N, alpha, X, incX, A, lda);
}

void her(const enum MATRIX_ORDER Order, const enum MATRIX_UPLO Uplo, const int N, const double
alpha, const std::complex<double> *X, const int incX, std::complex<double> *A, const int lda) {
	cblas_zher ((CBLAS_ORDER) Order, (CBLAS_UPLO) Uplo, N, alpha, X, incX, A, lda);
}

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc) {
	cblas_sgemm ((CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, (CBLAS_TRANSPOSE) TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
	cblas_dgemm ((CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, (CBLAS_TRANSPOSE) TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *B, const int ldb, const std::complex<float> beta, std::complex<float> *C, const int ldc) {
	cblas_cgemm ((CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, (CBLAS_TRANSPOSE) TransB, M, N, K, (const MKL_Complex8 *) &alpha, (const MKL_Complex8 *) A, lda, (const MKL_Complex8 *) B, ldb, (const MKL_Complex8 *) &beta, (MKL_Complex8 *) C, ldc);
}

void gemm(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const enum MATRIX_TRANSPOSE TransB, const int M, const int N, const int K, const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *B, const int ldb, const std::complex<double> beta, std::complex<double> *C, const int ldc) {
	cblas_zgemm ((CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, (CBLAS_TRANSPOSE) TransB, M, N, K, (const MKL_Complex16 *) &alpha, (const MKL_Complex16 *) A, lda, (const MKL_Complex16 *) B, ldb, (const MKL_Complex16 *) &beta, (MKL_Complex16 *) C, ldc);
}


void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const float alpha, const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY) {
	cblas_sgemv ( (CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY) {
	cblas_dgemv ( (CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const complex<float> alpha, const complex<float> *A, const int lda, const complex<float> *X, const int incX, const complex<float> beta, complex<float> *Y, const int incY) {
	cblas_cgemv ( (CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, M, N, &alpha, A, lda, X, incX, &beta, Y, incY);
}


void gemv(const enum MATRIX_ORDER Order, const enum MATRIX_TRANSPOSE TransA, const int M, const int N, const complex<double> alpha, const complex<double> *A, const int lda, const complex<double> *X, const int incX, const complex<double> beta, complex<double> *Y, const int incY) {
	cblas_zgemv ( (CBLAS_ORDER) Order, (CBLAS_TRANSPOSE) TransA, M, N, &alpha, A, lda, X, incX, &beta, Y, incY);
}

}
