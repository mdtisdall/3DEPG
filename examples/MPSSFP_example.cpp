#include "core/SpinStates.h"

#include "building_blocks/MPSSFPTR.h"

#include "BinaryFile.h"

#include <complex>

#include <iostream>

#include <cmath>

typedef float value_type;

void test() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates(10, 2.2e-3);

  MPSSFPTR<SpinStatesT, 2> mpssfpTR(
    value_type(10) * value_type(M_PI) / value_type(180),
    0,
    SpinStatesT::StateIndex({{1,0,0}}),
    4,
    1000,
    100,
    1e-4);

  { 
    SpinStatesT::cvalue_type excitedState = mpssfpTR(&spinStates); 

    std::cout << "MPSSFP TR signal: " << excitedState << std::endl;
    
    std::cout << "after MPSSFP TR:" << std::endl;
    for(int i = -1; i <= 1; i++) {
      SpinStatesT::StateIndex curIndex({{i,0,0}});
      SpinStatesT::SpinState curState = spinStates.getState(curIndex);
     
      std::cout << "\t(" << curIndex[0] << ", "
        << curIndex[1] << ", " << curIndex[2] << "): (" <<
        curState.fPlus << ", " <<
        curState.fMinus << ", " <<
        curState.z << ")" << std::endl;
    }
  }
}

int main() {
  test();

  return 0;
}

