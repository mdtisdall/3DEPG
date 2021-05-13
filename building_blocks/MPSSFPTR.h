#ifndef MPSSFPTR_h
#define MPSSFPTR_h

template <typename SpinStates, unsigned int NumPulses>
class MPSSFPTR {
  public:
    typedef typename SpinStates::value_type value_type;

    MPSSFPTR(
      const value_type flipAngle,
      const value_type flipPhaseStep,
      const typename SpinStates::StateIndex spoilingGrad,
      const value_type interPulseTime,
      const value_type t1,
      const value_type t2,
      const value_type threshold = value_type(0)
      ) :
      flipAngle(flipAngle),
      flipCurPhase(0),
      flipPhaseStep(flipPhaseStep),
      spoiler(spoilingGrad, interPulseTime, threshold),
      relaxation(interPulseTime, t1, t2)
    {}

    typename SpinStates::cvalue_type operator()(SpinStates *states) {
      for(unsigned int i = 0; i < NumPulses; i++) {
        { 
          typename SpinStates::Excitation excitation(flipAngle, flipCurPhase);
          flipCurPhase += flipPhaseStep;
          while(flipCurPhase > value_type(2.0 * M_PI) ) {
            flipCurPhase -= value_type(2.0 * M_PI);
          }
          while(flipCurPhase < value_type(-2.0 * M_PI) ) {
            flipCurPhase += value_type(2.0 * M_PI);
          }

          excitation(states);
        }
        spoiler(states);
        relaxation(states);
      } 

      typename SpinStates::cvalue_type ret =
        states->getState(typename SpinStates::StateIndex({{0,0,0}})).fPlus;
      spoiler(states);
      relaxation(states); 

      return ret;
    }


  protected:
    value_type flipAngle;
    value_type flipCurPhase;
    value_type flipPhaseStep;
    typename SpinStates::Spoiling spoiler;
    typename SpinStates::Relaxation relaxation;
};

#endif
