#include "SpinStates.h"

#include "FLASHTR.h"
#include "MPRAGE.h"

#include "ConstantGradSpoilingController.h"
#include "HexagonalGradSpoilingController.h"

#include "ConstantRFSpoilingController.h"
#include "ExponentialRFSpoilingController.h"
#include "ResettingExponentialRFSpoilingController.h"

#include "BinaryFile.h"

#include <complex>

#include <iostream>

#include <cmath>

typedef float value_type;

void test1() {
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
}

void test2() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates;

 typedef ConstantGradSpoilingController<SpinStatesT::StateIndex> 
   GradSpoilingController;

 typedef ConstantRFSpoilingController<SpinStatesT::value_type>
    RFSpoilingController;

  MPRAGE<
    SpinStatesT,
    GradSpoilingController,
    RFSpoilingController
  > mprage(
    value_type(7) * value_type(M_PI) / value_type(180),
    500,
    1200,
    176,
    256,
    6.12,
    1000,
    100,
    GradSpoilingController(),
    RFSpoilingController(),
    1e-4);

  
  { 
    std::vector<SpinStatesT::cvalue_type> measuredSignals;
    std::vector<value_type> rfPhases;
    mprage(&spinStates, &measuredSignals, &rfPhases); 
    
    std::cout << "MPRAGE signal length: " << measuredSignals.size()
      << std::endl;
    
    std::cout << "MPRAGE end signal: " << measuredSignals.back()
      << std::endl;
    
    std::cout << "MPRAGE end spin states: " << spinStates.size()
      << std::endl;

    BinaryFile<
      std::vector<SpinStatesT::cvalue_type>
      >::write(
      &measuredSignals, "test_constant_output.dat");

    BinaryFile<
      std::vector<value_type>
      >::write(
      &rfPhases, "test_constant_phases.dat");
  }
}

void test3() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates;


  typedef HexagonalGradSpoilingController<SpinStatesT::StateIndex> 
    GradSpoilingController;

  typedef ConstantRFSpoilingController<SpinStatesT::value_type>
    RFSpoilingController;

  MPRAGE<
    SpinStatesT,
    GradSpoilingController,
    RFSpoilingController
  > mprage(
    value_type(7) * value_type(M_PI) / value_type(180),
    500,
    1200,
    176,
    256,
    6.12,
    1000,
    100,
    GradSpoilingController(),
    RFSpoilingController(),
    1e-4);

  
  { 
    std::vector<SpinStatesT::cvalue_type> measuredSignals;
    std::vector<value_type> rfPhases;
    mprage(&spinStates, &measuredSignals, &rfPhases); 
    
    std::cout << "MPRAGE signal length: " << measuredSignals.size()
      << std::endl;
    
    std::cout << "MPRAGE end signal: " << measuredSignals.back()
      << std::endl;
    
    std::cout << "MPRAGE end spin states: " << spinStates.size()
      << std::endl;

    BinaryFile<
      std::vector<SpinStatesT::cvalue_type>
      >::write(
      &measuredSignals, "test_hexagonal_output.dat");

    BinaryFile<
      std::vector<value_type>
      >::write(
      &rfPhases, "test_hexagonal_phases.dat");
  }
}

void test4() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates;

 typedef ConstantGradSpoilingController<SpinStatesT::StateIndex> 
   GradSpoilingController;

 typedef ExponentialRFSpoilingController<SpinStatesT::value_type>
    RFSpoilingController;

  MPRAGE<
    SpinStatesT,
    GradSpoilingController,
    RFSpoilingController
  > mprage(
    value_type(7) * value_type(M_PI) / value_type(180),
    500,
    1200,
    176,
    256,
    6.12,
    1000,
    100,
    GradSpoilingController(),
    RFSpoilingController(value_type(50) * value_type(M_PI) / value_type(180)),
    1e-4);

  
  { 
    std::vector<SpinStatesT::cvalue_type> measuredSignals;
    std::vector<value_type> rfPhases;
    mprage(&spinStates, &measuredSignals, &rfPhases); 
    
    std::cout << "MPRAGE signal length: " << measuredSignals.size()
      << std::endl;
    
    std::cout << "MPRAGE end signal: " << measuredSignals.back()
      << std::endl;
    
    std::cout << "MPRAGE end spin states: " << spinStates.size()
      << std::endl;

    BinaryFile<
      std::vector<SpinStatesT::cvalue_type>
      >::write(
      &measuredSignals, "test_constant_rfspoil_output.dat");

    BinaryFile<
      std::vector<value_type>
      >::write(
      &rfPhases, "test_constant_rfspoil_phases.dat");
  }
}

void test5() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates;


  typedef HexagonalGradSpoilingController<SpinStatesT::StateIndex> 
    GradSpoilingController;

  typedef ExponentialRFSpoilingController<SpinStatesT::value_type>
    RFSpoilingController;

  MPRAGE<
    SpinStatesT,
    GradSpoilingController,
    RFSpoilingController
  > mprage(
    value_type(7) * value_type(M_PI) / value_type(180),
    500,
    1200,
    176,
    256,
    6.12,
    1000,
    100,
    GradSpoilingController(),
    RFSpoilingController(value_type(50) * value_type(M_PI) / value_type(180)),
    1e-4);

  
  { 
    std::vector<SpinStatesT::cvalue_type> measuredSignals;
    std::vector<value_type> rfPhases;
    mprage(&spinStates, &measuredSignals, &rfPhases); 
    
    std::cout << "MPRAGE signal length: " << measuredSignals.size()
      << std::endl;
    
    std::cout << "MPRAGE end signal: " << measuredSignals.back()
      << std::endl;
    
    std::cout << "MPRAGE end spin states: " << spinStates.size()
      << std::endl;

    BinaryFile<
      std::vector<SpinStatesT::cvalue_type>
      >::write(
      &measuredSignals, "test_hexagonal_rfspoil_output.dat");

    BinaryFile<
      std::vector<value_type>
      >::write(
      &rfPhases, "test_hexagonal_rfspoil_phases.dat");
  }
}

void test6() {
  typedef SpinStates<long, value_type> SpinStatesT;
  SpinStatesT spinStates;


  typedef HexagonalGradSpoilingController<SpinStatesT::StateIndex> 
    GradSpoilingController;

  typedef ResettingExponentialRFSpoilingController<SpinStatesT::value_type>
    RFSpoilingController;

  MPRAGE<
    SpinStatesT,
    GradSpoilingController,
    RFSpoilingController
  > mprage(
    value_type(7) * value_type(M_PI) / value_type(180),
    500,
    1200,
    176,
    256,
    6.12,
    1000,
    100,
    GradSpoilingController(),
    RFSpoilingController(
      value_type(50) * value_type(M_PI) / value_type(180),
      176),
    1e-4);

  
  { 
    std::vector<SpinStatesT::cvalue_type> measuredSignals;
    std::vector<value_type> rfPhases;
    mprage(&spinStates, &measuredSignals, &rfPhases); 
    
    std::cout << "MPRAGE signal length: " << measuredSignals.size()
      << std::endl;
    
    std::cout << "MPRAGE end signal: " << measuredSignals.back()
      << std::endl;
    
    std::cout << "MPRAGE end spin states: " << spinStates.size()
      << std::endl;

    BinaryFile<
      std::vector<SpinStatesT::cvalue_type>
      >::write(
      &measuredSignals, "test_hexagonal_resetrfspoil_output.dat");

    BinaryFile<
      std::vector<value_type>
      >::write(
      &rfPhases, "test_hexagonal_resetrfspoil_phases.dat");
  }
}

int main() {

  test1();
  test2();
  test3();
  test4();
  test5();
  test6();

  return 0;
}
