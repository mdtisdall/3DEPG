#ifndef FLASHTR_h
#define FLASHTR_h

template <typename SpinStates>
class FLASHTR {
  public:
    typedef typename SpinStates::value_type value_type;

    FLASHTR(
      const value_type flipAngle,
      const value_type flipPhase,
      const typename SpinStates::StateIndex spoilingGrad,
      const value_type tr,
      const value_type t1,
      const value_type t2) :
      excitation(flipAngle, flipPhase),
      spoiling(spoilingGrad),
      relaxation(tr, t1, t2)
    {}

    void operator()(SpinStates *states) {
      excitation(states); 
      relaxation(states); 
      spoiling(states); 
    }


  protected:
    typename SpinStates::Excitation excitation;
    typename SpinStates::Spoiling spoiling;
    typename SpinStates::Relaxation relaxation;
};

#endif
