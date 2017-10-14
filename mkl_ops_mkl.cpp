#include "mkl_ops.h"

#include "mkl.h"

namespace MKL {

void vMul(const int n, const float *a, const float *b, float *y)
{
  vsMul(n, a, b, y);
}

void vMul(const int n, const double *a, const double *b, double *y)
{
  vdMul(n, a, b, y);
}

void vMul(const int n,
    const std::complex<float> *a,
    const std::complex<float> *b,
    std::complex<float> *y)
{
  vcMul(n, (MKL_Complex8*)a, (MKL_Complex8*)b, (MKL_Complex8*)y); 
}

void vMul(const int n,
    const std::complex<double> *a,
    const std::complex<double> *b,
    std::complex<double> *y)
{
  vzMul(n, (MKL_Complex16*)a, (MKL_Complex16*)b, (MKL_Complex16*)y); 
}

void vSqrt(const int n, const float *a, float *y)
{
  vsSqrt(n, a, y);
}

void vSqrt(const int n, const double *a, double *y)
{
  vdSqrt(n, a, y);
}

void vSqrt(const int n, const std::complex<float> *a, std::complex<float> *y)
{
  vcSqrt(n, (MKL_Complex8*)a, (MKL_Complex8*)y);
}

void vSqrt(const int n, const std::complex<double> *a, std::complex<double> *y)
{
  vzSqrt(n, (MKL_Complex16*)a, (MKL_Complex16*)y);
}


void vExp(const int n, const float *a, float *y)
{
  vsExp(n, a, y);
}

void vExp(const int n, const double *a, double *y)
{
  vdExp(n, a, y);
}

void vExp(const int n, const std::complex<float> *a, std::complex<float> *y)
{
  vcExp(n, (MKL_Complex8*)a, (MKL_Complex8*)y);
}

void vExp(const int n, const std::complex<double> *a, std::complex<double> *y)
{
  vzExp(n, (MKL_Complex16*)a, (MKL_Complex16*)y);
}

}
