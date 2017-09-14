#include "SpinStates.h"

#include <complex>

#include <iostream>

int main() {
  typedef SpinStates<long, float> SpinStatesT;
  SpinStatesT spinStates;

  {
    SpinStatesT::StateIndex origin = {0,0,0};
    SpinStatesT::SpinState originState = spinStates.getState(origin);

    std::cout << "before excitation: (" <<
      originState.fPlus << ", " <<
      originState.fMinus << ", " <<
      originState.z << ")" << std::endl;
  }

  SpinStatesT::Excitation excitation(1.0, 0.0);
  excitation.excite(&spinStates);
  
  {  
    SpinStatesT::StateIndex origin = {0,0,0};
    SpinStatesT::SpinState originState = spinStates.getState(origin);

    std::cout << "after excitation: (" <<
      originState.fPlus << ", " <<
      originState.fMinus << ", " <<
      originState.z << ")" << std::endl; 
  }

  SpinStatesT::Relaxation relaxation(100, 1000, 100);
  relaxation.relax(&spinStates);
  
  {  
    SpinStatesT::StateIndex origin = {0,0,0};
    SpinStatesT::SpinState originState = spinStates.getState(origin);
    
    std::cout << "after relaxation: (" <<
      originState.fPlus << ", " <<
      originState.fMinus << ", " <<
      originState.z << ")" << std::endl;
  }
  
  { 
    SpinStatesT::Spoiling spoiling(SpinStatesT::StateIndex({1,0,0}));
    spoiling.spoil(&spinStates); 

    std::cout << "after spoiling:" << std::endl;
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
  
  excitation.excite(&spinStates);
  relaxation.relax(&spinStates);
  
  { 
    SpinStatesT::Spoiling spoiling(SpinStatesT::StateIndex({0,1,0}));
    spoiling.spoil(&spinStates); 

    std::cout << "after spoiling:" << std::endl;
    for(int z = -1; z <= 1; z++) {
      for(int y = -1; y <= 1; y++) {
        SpinStatesT::StateIndex curIndex = {z,y,0};
        SpinStatesT::SpinState curState = spinStates.getState(curIndex);
     
        std::cout << "\t(" << curIndex[0] << ", "
          << curIndex[1] << ", " << curIndex[2] << "): (" <<
          curState.fPlus << ", " <<
          curState.fMinus << ", " <<
          curState.z << ")" << std::endl;
      }
    }
  }
  
  { 
    SpinStatesT::Spoiling spoiling(SpinStatesT::StateIndex({0,0,1}), 0.05);
    spoiling.spoil(&spinStates); 

    std::cout << "after spoiling:" << std::endl;
    for(int z = -1; z <= 1; z++) {
      for(int y = -1; y <= 1; y++) {
        for(int x = -1; x <= 1; x++) {
          SpinStatesT::StateIndex curIndex = {z,y,x};
          SpinStatesT::SpinState curState = spinStates.getState(curIndex);
     
          std::cout << "\t(" << curIndex[0] << ", "
            << curIndex[1] << ", " << curIndex[2] << "): (" <<
            curState.fPlus << ", " <<
            curState.fMinus << ", " <<
            curState.z << ")" << std::endl;
        }
      }
    }
  }


  return 0;
}
