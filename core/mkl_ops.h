#ifndef mkl_ops_h
#define mkl_ops_h

#include<complex>

namespace MKL {

void vMul(const int n, const float *a, const float *b, float *y); 
void vMul(const int n, const double *a, const double *b, double *y); 
void vMul(const int n,
    const std::complex<float> *a,
    const std::complex<float> *b,
    std::complex<float> *y); 
void vMul(const int n,
    const std::complex<double> *a,
    const std::complex<double> *b,
    std::complex<double> *y); 

void vSqrt(const int n, const float *a, float *y); 
void vSqrt(const int n, const double *a, double *y); 
void vSqrt(const int n, const std::complex<float> *a, std::complex<float> *y); 
void vSqrt(const int n, const std::complex<double> *a, std::complex<double> *y); 

void vExp(const int n, const float *a, float *y); 
void vExp(const int n, const double *a, double *y); 
void vExp(const int n, const std::complex<float> *a, std::complex<float> *y); 
void vExp(const int n, const std::complex<double> *a, std::complex<double> *y); 

}

#endif
