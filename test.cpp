#include "SpinStates.h"

#include "FLASHTR.h"

#include <complex>

#include <iostream>

#include <cmath>

typedef float value_type;

int main() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates;

  FLASHTR<SpinStatesT> flashTR(
    value_type(7) * value_type(M_PI) / value_type(180),
    0,
    SpinStatesT::StateIndex({1,0,0}),
    10,
    1000,
    100);

  
  { 
    SpinStatesT::cvalue_type excitedState = flashTR(&spinStates); 

    std::cout << "FLASH TR signal: " << excitedState << std::endl;
    
    std::cout << "after FLASH TR:" << std::endl;
    for(int i = -1; i <= 1; i++) {
      SpinStatesT::StateIndex curIndex = {i,0,0};
      SpinStatesT::SpinState curState = spinStates.getState(curIndex);
     
      std::cout << "\t(" << curIndex[0] << ", "
        << curIndex[1] << ", " << curIndex[2] << "): (" <<
        curState.fPlus << ", " <<
        curState.fMinus << ", " <<
        curState.z << ")" << std::endl;
    }
  }

  return 0;
}
