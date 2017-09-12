#include "SpinStates.h"

#include <complex>

#include <iostream>

int main() {
  typedef SpinStates<long, float> SpinStatesT;
  SpinStatesT spinStates;

  SpinStatesT::StateIndex origin(0,0,0);
  SpinStatesT::SpinState originState = spinStates.getState(origin);

  std::cout << "before excitation: (" <<
    originState.fPlus << ", " <<
    originState.fMinus << ", " <<
    originState.z << ")" << std::endl; 

  SpinStatesT::Excitation excitation(1.0, 0.0);
  excitation.excite(&spinStates);
  
  originState = spinStates.getState(origin);

  std::cout << "after excitation: (" <<
    originState.fPlus << ", " <<
    originState.fMinus << ", " <<
    originState.z << ")" << std::endl; 
  return 0;
}
