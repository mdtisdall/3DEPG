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

  std::array<value_type, 3> flipAngles;

  for(unsigned int fa = 1; fa < 180; fa++) {
    flipAngles.fill(value_type(fa) * value_type(M_PI) / value_type(180));
    //flipAngles[0] = value_type(10) * value_type(M_PI) / value_type(180);
    //flipAngles[1] = value_type(30) * value_type(M_PI) / value_type(180);
    //flipAngles[2] = value_type(10) * value_type(M_PI) / value_type(180);

    MPSSFPTR<SpinStatesT, 3> mpssfpTR(
      flipAngles,
      //0,
      value_type(180) * value_type(M_PI) / value_type(180),
      SpinStatesT::StateIndex({{1,0,0}}),
      4,
      650,
      80,
      value_type(0) * value_type(M_PI) / value_type(180),
      1e-4);

    {
      for(unsigned int i = 0; i < 1000; i++) {
        mpssfpTR(&spinStates);
      }

      SpinStatesT::cvalue_type excitedState = mpssfpTR(&spinStates); 

      std::cout << ", " << std::abs(excitedState) << std::endl;
      
      //std::cout << "after MPSSFP TR:" << std::endl;
      //for(int i = -1; i <= 1; i++) {
      //  SpinStatesT::StateIndex curIndex({{i,0,0}});
      //  SpinStatesT::SpinState curState = spinStates.getState(curIndex);
      // 
      //  std::cout << "\t(" << curIndex[0] << ", "
      //    << curIndex[1] << ", " << curIndex[2] << "): (" <<
      //    curState.fPlus << ", " <<
      //    curState.fMinus << ", " <<
      //    curState.z << ")" << std::endl;
      //}
    }
  }
}

int main() {
  test();

  return 0;
}

