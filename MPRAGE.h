#ifndef MPRAGE_h
#define MPRAGE_h

#include "FLASHTR.h"

template <
  typename SpinStates,
  typename GradSpoilingController,
  typename RFSpoilingController>
class MPRAGE {
  public:
    typedef typename SpinStates::value_type value_type;

    MPRAGE(
      const value_type flipAngle,
      const value_type td1,
      const value_type td2,
      const unsigned int innerSteps,
      const unsigned int outerSteps,
      const value_type flashTR,
      const value_type t1,
      const value_type t2,
      const GradSpoilingController gradSpoilingController,
      const RFSpoilingController rfSpoilingController,
      const value_type threshold = value_type(0)
    ) :
    flipAngle(flipAngle),
    innerSteps(innerSteps),
    outerSteps(outerSteps),
    flashTR(flashTR),
    t1(t1),
    t2(t2),
    inversion(value_type(M_PI), value_type(0)),
    td1Relaxation(td1, t1, t2),
    td2Relaxation(td2, t1, t2),
    gradSpoilingController(gradSpoilingController),
    rfSpoilingController(rfSpoilingController),
    threshold(threshold)
  {}

    void operator()(
      SpinStates *states,
      std::vector<typename SpinStates::cvalue_type> *output) {
      for(unsigned int outerStep = 0; outerStep < outerSteps; outerStep++) {
        inversion(states);
        td1Relaxation(states);

        for(unsigned int innerStep = 0; innerStep < innerSteps; innerStep++) {
          output->push_back(
            FLASHTR<SpinStates>(
              flipAngle, rfSpoilingController.nextPhase(),
              gradSpoilingController.nextGrad(),
              flashTR, t1, t2, threshold)(states)
          );
        }
        
        td2Relaxation(states);
      }
    }
  protected:
    const value_type flipAngle;
    const unsigned int innerSteps;
    const unsigned int outerSteps;
    const value_type flashTR;
    const value_type t1;
    const value_type t2;
    typename SpinStates::Excitation inversion;
    typename SpinStates::Relaxation td1Relaxation;
    typename SpinStates::Relaxation td2Relaxation;
    GradSpoilingController gradSpoilingController;
    RFSpoilingController rfSpoilingController;
    const value_type threshold;
};

#endif
