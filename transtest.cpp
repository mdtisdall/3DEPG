
#include "mkl.h"
#include <iostream>
#include <chrono>

using namespace std;

void form1(
  float *states,
  float *trans,
  float *output) {
  
  auto start = chrono::steady_clock::now();

  for(unsigned int i = 0; i < 100; i ++) {
    cblas_sgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        3, 256 * 256 * 256, 3,
        1.0, trans, 3, states, 3, 0, output, 3);
  }
  
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

void form2(
  float *states,
  float *trans,
  float *output) {
  
  auto start = chrono::steady_clock::now();

  for(unsigned int i = 0; i < 100; i ++) {
    cblas_sgemm(
      CblasColMajor, CblasNoTrans, CblasTrans,
      3, 256 * 256 * 256, 3,
      1.0, trans, 3, states, 256 * 256 * 256, 0, output, 3);
  }
  
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

void form3(
  float *states,
  float *trans,
  float *output) {
  
  auto start = chrono::steady_clock::now();

  for(unsigned int i = 0; i < 100; i ++) {
    cblas_sgemm(
      CblasColMajor, CblasNoTrans, CblasNoTrans,
      256 * 256 * 256, 3, 3,
      1.0, states, 256 * 256 * 256, trans, 3, 0, output, 256 * 256 * 256);
  }
  
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

int main() {
  float *states = (float*) malloc(sizeof(float) * 3 * 256 * 256 * 256);
  float trans [3 * 3];
  float *output = (float*) malloc(sizeof(float) * 3 * 256 * 256 * 256);
  
   
  form1(states, trans, output);
  form2(states, trans, output);
  form3(states, trans, output);
  

  return 0;
}
